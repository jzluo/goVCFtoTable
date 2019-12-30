package main

import (
	"bufio"
	"compress/gzip"
	"fmt"
	"os"
	"strings"
)

type Variant struct {
	Chr, Pos int
	Ref, Alt string
	Calls    Call
}

type Call struct {
	PT_ID, Var_ID, RNC string
	GT, GQ, Ref, Alt   int
}

func ReadVCFgz(filename string) (*bufio.Reader, []string) {
	f, err := os.Open(filename)
	if err != nil {
		fmt.Println("File not found")
		os.Exit(0)
	}
	defer f.Close()

	gz, err := gzip.NewReader(f)
	if err != nil {
		fmt.Println("Error reading file")
		os.Exit(1)
	}
	defer gz.Close()

	reader := bufio.NewReader(gz)
	var sampleNames []string

	for {
		line, err := reader.ReadString('\n')
		if err != nil {
			fmt.Println("Error parsing VCF file")
			fmt.Println(err)
			os.Exit(3)
		}

		if strings.HasPrefix(line, "#CHROM") {
			sampleNames = ParseSampleNames(line)
			break
		}
	}

	return reader, sampleNames
}

func ParseSampleNames(line string) []string {
	var sampleNames []string
	splitLine := strings.Split(line, "\t")
	sampleNames = splitLine[9:]
	for i := range sampleNames {
		sampleNames[i] = strings.Split(sampleNames[i], "_")[1]
	}
	return sampleNames
}

// GT:DP:AD:SBPV:GQ:PL:FT:RNC

func main() {
	_, sampleNames := ReadVCFgz("PKD1_IDT_F145K.vcf.gz")
	fmt.Println(sampleNames)
}
