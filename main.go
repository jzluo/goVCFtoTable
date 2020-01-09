package main

import (
	"bufio"
	"bytes"
	"flag"
	"fmt"
	"github.com/biogo/hts/bgzf"
	"io"
	"log"
	"os"
	"runtime"
	"runtime/pprof"
	"time"
)

func CheckError(err error) bool {
	if err != nil {
		if err == io.EOF {
			return true
		}
	} else {
		fmt.Println(err)
		os.Exit(0)
	}
	return false
}

//func ReadVCFgz(filename string) (*bufio.Reader, *os.File, *gzip.Reader, []string) {
func ReadVCFgz(filename string) (*bufio.Reader, *os.File, *bgzf.Reader, [][]byte, int) {
	f, err := os.Open(filename)
	if err != nil {
		fmt.Println("File not found")
		os.Exit(0)
	}

	gz, err := bgzf.NewReader(f, 8)
	if err != nil {
		fmt.Println("Error reading file")
		os.Exit(1)
	}

	reader := bufio.NewReader(gz)

	for {
		line, err := reader.ReadBytes('\n')
		if err != nil {
			fmt.Println("Error parsing VCF file")
			fmt.Println(err)
			os.Exit(3)
		}

		if bytes.HasPrefix(line, []byte("#CHROM")) {
			sampleNames, numSamples := ParseSampleNames(line)
			return reader, f, gz, sampleNames, numSamples
		}
	}
}

func ParseSampleNames(line []byte) ([][]byte, int) {
	sampleNames := make([][]byte, 0)
	//splitLine := bytes.Split(line, []byte("\t"))
	var start int
	_, start = GetNextCol(line, start, 9, '\t')

	for start > 0 { //GetNextCol() returns 0 when there are no more delimiters
		var sample []byte
		sample, start = GetNextCol(line, start, 0, '\t')
		sampleNames = append(sampleNames, sample)
	}

	return sampleNames, len(sampleNames)
}

func GetNextCol(row []byte, start int, skip int, delim byte) ([]byte, int) {
	for i := 0; i < skip; i++ {
		start = bytes.IndexByte(row[start:], delim) + start + 1
	}
	nextTabPos := bytes.IndexByte(row[start:], delim)
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

func VCFtoTable(reader *bufio.Reader, writer *bufio.Writer, sampleNames [][]byte, numSamples int) {
	for {
		line, err := reader.ReadBytes('\n')
		if err == io.EOF {
			break
		} else if err != nil {
			fmt.Println(err)
			os.Exit(1)
		}

		ProcessVariant(line, writer, sampleNames, numSamples)
	}
}

func ProcessVariant(row []byte, writer *bufio.Writer, sampleNames [][]byte, numSamples int) {
	//var chr, pos, ref, alt []byte
	//var start int
	//chr, start = GetNextCol(row, 0, 0, '\t')
	//pos, start = GetNextCol(row, start, 0, '\t')
	//ref, start = GetNextCol(row, start, 1, '\t')
	//alt, start = GetNextCol(row, start, 0, '\t')
	_, start := GetNextCol(row, 0, 9, '\t')
	ProcessCalls(row[start:], writer, sampleNames, numSamples)
}

func ProcessCalls(calls []byte, writer *bufio.Writer, sampleNames [][]byte, numSamples int) {
	var start int
	het := []byte("0/1")
	hetcmp := []byte("compound")
	hom := []byte("1/1")
	com := []byte(",")
	newline := []byte("\n")
	for i := 0; i < numSamples; i++ {
		var call, flag, gt []byte
		call, start = GetNextCol(calls, start, 0, '\t')

		if call[0] == '.' || call[2] == '.' {
			continue
		}

		if call[0] == 48 && call[2] == 49 {
			gt = het
		} else if call[0] == 49 && call[2] == 49 {
			gt = hom
		} else if call[0] == 49 && call[2] == 48 {
			gt = het
			flag = hetcmp
		} else if call[0] == 48 && call[2] == 48 {
			continue
		} else {
			fmt.Println(string(call[:3]))
			panic("wtfff")
		}
		dp, ad_ref, ad_alt, gq, ft := GetFormatFields(call)

		writer.Write(bytes.Join([][]byte{
			sampleNames[i],
			gt,
			gq,
			ad_ref,
			ad_alt,
			dp,
			ft,
			flag,
			newline,
		}, com))
	}
}

var cpuprofile = flag.String("cpuprofile", "", "write cpu profile to `file`")
var memprofile = flag.String("memprofile", "", "write memory profile to `file`")

func main() {
	flag.Parse()
	if *cpuprofile != "" {
		f, err := os.Create(*cpuprofile)
		if err != nil {
			log.Fatal("could not create CPU profile: ", err)
		}
		defer f.Close()
		if err := pprof.StartCPUProfile(f); err != nil {
			log.Fatal("could not start CPU profile: ", err)
		}
		defer pprof.StopCPUProfile()
	}
	now := time.Now()
	fmt.Println(now)
	reader, f, gz, sampleNames, numSamples := ReadVCFgz("PKD1_IDT_F145K.vcf.gz")
	file, err := os.Create("testout_fancy.csv")
	if err != nil {
		fmt.Println(err)
		os.Exit(5)
	}

	writer := bufio.NewWriter(file)
	defer writer.Flush()
	defer f.Close()
	defer gz.Close()

	VCFtoTable(reader, writer, sampleNames, numSamples)

	fmt.Println(time.Since(now))

	if *memprofile != "" {
		f, err := os.Create(*memprofile)
		if err != nil {
			log.Fatal("could not create memory profile: ", err)
		}
		defer f.Close()
		runtime.GC() // get up-to-date statistics
		if err := pprof.WriteHeapProfile(f); err != nil {
			log.Fatal("could not write memory profile: ", err)
		}
	}

}
