package main

import (
	"bufio"
	"io"
	"log"
	"os"

	arg "github.com/alexflint/go-arg"
	"github.com/biogo/hts/bam"
	"github.com/biogo/hts/sam"
	"github.com/brentp/xopen"
)

// VERSION defines the program version.
const VERSION = "0.2"

// Opts is the struct with the options that the program accepts.
// Opts encapsulates common command line options.
type Opts struct {
	Input  string `arg:"positional,required" help:"file (- for STDIN)"`
	Tag    string `arg:"required" help:"tag to be added with count"`
	Sorted bool   `arg:"required" help:"file is sorted by QNAME (use for forward compatibility)"`
}

// Version returns the program name and version.
func (Opts) Version() string { return "sam-count-secondary " + VERSION }

// Description returns an extended description of the program.
func (Opts) Description() string {
	return "Counts and adds a tag with the number of alignments for each read. Currently it supports only BAM files sorted by QNAME. Rules are not well defined as to how supplementary/chimeric alignments should be handled. Currently chimeric alignments are handled like any other alignment and are counted as secondary"
}

func main() {
	var opts Opts
	opts.Tag = "X0"
	arg.MustParse(&opts)

	// Open file for reading.
	var fh io.Reader
	var err error
	if opts.Input == "-" {
		if !xopen.IsStdin() {
			log.Fatal("cannot read from stdin")
		}
		fh = os.Stdin
	} else {
		fh, err = os.Open(opts.Input)
		if err != nil {
			log.Fatalf("cannot open file: %v", err)
		}
	}

	// Open SAM/BAM reader.
	br, err := bam.NewReader(fh, 2)
	if err != nil {
		log.Fatalf("cannot create sam reader: %v", err)
	}

	// Close all readers at the end.
	defer func() {
		if err = br.Close(); err != nil {
			log.Fatalf("cannot close sam reader: %v", err)
		}
	}()

	// Create a SAM writer that writes to STDOUT.
	stdout := bufio.NewWriter(os.Stdout)
	defer stdout.Flush()
	w, err := sam.NewWriter(stdout, br.Header(), sam.FlagDecimal)
	if err != nil {
		log.Fatalf("write of header failed: %v", err)
	}

	// Loop on the SAM/BAM file.
	samBlk := make([]*sam.Record, 0)
	currName := ""
	warnedOnce := false
	for {
		r, err := br.Read()
		if err != nil {
			if err == io.EOF {
				break
			}
			log.Fatalf("error reading SAM: %v", err)
		}
		// Warn for supplementary/chimeric alignments because it is not well
		// defined how these should be handled.
		if _, ok := r.Tag([]byte("SA")); ok && !warnedOnce {
			log.Print("warning: not optimized to handle supplementary/chimeric alignments; consider filtering out records with SA tag")
			warnedOnce = true
		}

		if currName == "" || r.Name == currName {
			samBlk = append(samBlk, r)
			currName = r.Name
			continue
		}

		// Process previous block.
		processBlk(samBlk, w, opts.Tag)

		// Create new block
		samBlk = nil
		samBlk = append(samBlk, r)
		currName = r.Name
	}
	processBlk(samBlk, w, opts.Tag)
}

func processBlk(samBlk []*sam.Record, w *sam.Writer, tag string) {
	cnt := len(samBlk)
	for _, pr := range samBlk {
		aux, err := sam.NewAux(sam.NewTag(tag), cnt)
		if err != nil {
			log.Fatalf("failed to create aux tag for %s: %v ", pr.Name, err)
		}
		pr.AuxFields = append(pr.AuxFields, aux)

		if err := w.Write(pr); err != nil {
			log.Fatalf("write failed: %v for %s", err, pr.Name)
		}
	}
}
