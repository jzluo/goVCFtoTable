package main

import (
	"bufio"
	"bytes"
	"strconv"
	//"flag"
	"fmt"
	"github.com/biogo/hts/bgzf"
	"io"
	//"log"
	"os"
	//"runtime"
	//"runtime/pprof"
	"time"
)

//var cpuprofile = flag.String("cpuprofile", "", "write cpu profile to `file`")
//var memprofile = flag.String("memprofile", "", "write memory profile to `file`")

func main() {
	//flag.Parse()
	//if *cpuprofile != "" {
	//	f, err := os.Create(*cpuprofile)
	//	if err != nil {
	//		log.Fatal("could not create CPU profile: ", err)
	//	}
	//	defer f.Close()
	//	if err := pprof.StartCPUProfile(f); err != nil {
	//		log.Fatal("could not start CPU profile: ", err)
	//	}
	//	defer pprof.StopCPUProfile()
	//}
	now := time.Now()
	fileName := os.Args[1]
	outFile := os.Args[2]
	fmt.Println(now, "Doing", fileName)
	threads, err := strconv.Atoi(os.Args[3])
	if err != nil {
		fmt.Println("invalid number of threads")
		os.Exit(9)
	}

	VCFtoTable(fileName, outFile, threads)

	fmt.Println("Finished", fileName, "in", time.Since(now))

	//if *memprofile != "" {
	//	f, err := os.Create(*memprofile)
	//	if err != nil {
	//		log.Fatal("could not create memory profile: ", err)
	//	}
	//	defer f.Close()
	//	runtime.GC() // get up-to-date statistics
	//	if err := pprof.WriteHeapProfile(f); err != nil {
	//		log.Fatal("could not write memory profile: ", err)
	//	}
	//}
}

func ReadVCFgz(filename string, threads int) (*bufio.Reader, *os.File, *bgzf.Reader, [][]byte, int) {
	f, err := os.Open(filename)
	if err != nil {
		fmt.Println("File not found")
		os.Exit(0)
	}

	gz, err := bgzf.NewReader(f, threads)
	if err != nil {
		fmt.Println("Error decompressing file")
		os.Exit(1)
	}

	reader := bufio.NewReader(gz)

	for {
		// Read until the last header line and return samples
		line, err := reader.ReadBytes('\n')
		if err != nil {
			fmt.Println("Error parsing VCF file")
			fmt.Println(err)
			os.Exit(3)
		}

		if bytes.HasPrefix(line, []byte("#CHROM")) {
			line = line[:len(line)-1]
			sampleNames, numSamples := ParseSampleNames(line)
			return reader, f, gz, sampleNames, numSamples
		}
	}
}

func VCFtoTable(fileName, outFile string, threads int) {
	reader, f, gz, sampleNames, numSamples := ReadVCFgz(fileName, threads)
	file, err := os.Create(outFile)
	if err != nil {
		fmt.Println(err)
		os.Exit(5)
	}

	writer := bufio.NewWriter(file)
	defer writer.Flush()
	defer gz.Close()
	defer f.Close()

	//writer.WriteString("Var_ID\tPT_ID\tGT\tGQ\tAD_REF\tAD_ALT\tDP\tFT\n") // write header

	for {
		line, err := reader.ReadBytes('\n')
		if err == io.EOF {
			break
		} else if err != nil {
			fmt.Println(err)
			os.Exit(1)
		}
		line = line[:len(line)-1] // strip \n
		ProcessVariant(line, writer, sampleNames, numSamples)
	}
}

func ProcessVariant(row []byte, writer *bufio.Writer, sampleNames [][]byte, numSamples int) {
	var varID []byte
	var start int
	//chr, start = GetNextCol(row, 0, 0, '\t')
	//pos, start = GetNextCol(row, start, 0, '\t')
	varID, start = GetNextCol(row, start, 2, '\t')
	//ref, start = GetNextCol(row, start, 0, '\t')
	//alt, start = GetNextCol(row, start, 0, '\t')
	_, start = GetNextCol(row, 0, 9, '\t')
	ProcessCalls(row[start:], writer, sampleNames, numSamples, varID)
}

// ProcessCalls loops through every sample's call and skips WT or no-calls (including monoallelic sites from GLnexus atm)
//
func ProcessCalls(calls []byte, writer *bufio.Writer, sampleNames [][]byte, numSamples int, varID []byte) {
	var start int
	het := []byte{49}
	cmp := []byte("CMPND")
	cmp2 := []byte(";CMPND")
	hom := []byte{50}
	tab := []byte("\t")
	newline := []byte("\n")
	zero := []byte{48}

	for i := 0; i < numSamples; i++ {
		var call, gt []byte
		compound := false
		call, start = GetNextCol(calls, start, 0, '\t')

		if call[0] == 48 && call[2] == 48 { // 0/0
			continue
		} else if call[0] == '.' || call[2] == '.' {
			continue
		} else if call[0] == 48 && call[2] == 49 { // 0/1
			gt = het
		} else if call[0] == 49 && call[2] == 49 { // 1/1
			gt = hom
		} else if call[0] == 49 && call[2] == 48 { // 1/0
			gt = het
			compound = true
		} else {
			fmt.Println(string(call[:3]))
			fmt.Println(string(varID))
			panic("wtf is this")
		}

		dp, ad_ref, ad_alt, gq, ft := GetFormatFields(call)
		if ad_ref[0] == '.' {
			ad_ref = zero
		}
		if compound {
			// allocate new slice to avoid overwriting of the original backing array (calls)
			ft = ft[:len(ft):len(ft)]
			if len(ft) == 0 {
				ft = cmp
			} else if len(ft) > 0 {
				ft = append(ft, cmp2...)
			}
		}

		writer.Write(bytes.Join([][]byte{
			varID,
			sampleNames[i],
			gt,
			gq,
			ad_ref,
			ad_alt,
			dp,
			ft,
		}, tab))
		writer.Write(newline)
	}
}

func ParseSampleNames(line []byte) ([][]byte, int) {
	sampleNames := make([][]byte, 0)
	var start int
	_, start = GetNextCol(line, start, 9, '\t')

	for start > 0 { //GetNextCol() returns 0 when there are no more delimiters
		var sample []byte
		sample, start = GetNextCol(line, start, 0, '\t')
		sampleNames = append(sampleNames, sample)
	}

	return sampleNames, len(sampleNames)
}

// GetNextCol returns the index of the next delimiter and the subslice from a starting index to the delimiter.
// The skip parameter allows for the skipping of a number of "columns"
func GetNextCol(row []byte, start int, skip int, delim byte) ([]byte, int) {
	for i := 0; i < skip; i++ {
		start = bytes.IndexByte(row[start:], delim) + start + 1
	}
	nextTabPos := bytes.IndexByte(row[start:], delim)
	// if end of line reached, return the remainder of the byte slice
	if nextTabPos == -1 {
		return row[start:], 0
	}
	return row[start : start+nextTabPos], start + nextTabPos + 1
}

func GetReads(ad []byte) ([]byte, []byte) {
	split := bytes.IndexByte(ad, ',')
	return ad[:split], ad[split+1:]
}

func GetFormatFields(call []byte) ([]byte, []byte, []byte, []byte, []byte) {
	// GT:DP:AD:SBPV:GQ:PL:FT:RNC
	var dp, ad, ad_ref, ad_alt, gq, ft []byte
	var start int
	dp, start = GetNextCol(call, 4, 0, ':')
	ad, start = GetNextCol(call, start, 0, ':')
	ad_ref, ad_alt = GetReads(ad)
	gq, start = GetNextCol(call, start, 1, ':')
	ft, start = GetNextCol(call, start, 1, ':')
	if ft[0] == 46 { // if ft=.
		ft = []byte("")
	}

	return dp, ad_ref, ad_alt, gq, ft
}
