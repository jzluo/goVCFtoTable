package main

import (
	"bufio"
	"compress/gzip"
	"encoding/csv"
	"flag"
	"fmt"
	"io"
	"log"
	"os"
	"runtime"
	"runtime/pprof"
	"strings"
	"time"
)

type Variant struct {
	Chr, Pos         string
	Ref, Alt, Var_ID string
	Calls            []Call
}

type Call struct {
	PT_ID, RNC             string
	GT, GQ, AD_Ref, AD_Alt string
}

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

func ReadVCFgz(filename string) (*bufio.Reader, *os.File, *gzip.Reader, []string) {
	f, err := os.Open(filename)
	if err != nil {
		fmt.Println("File not found")
		os.Exit(0)
	}
	//defer f.Close()

	gz, err := gzip.NewReader(f)
	if err != nil {
		fmt.Println("Error reading file")
		os.Exit(1)
	}
	//defer gz.Close()

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

	return reader, f, gz, sampleNames
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
	reader, f, gz, sampleNames := ReadVCFgz("PKD1_IDT_F145K.vcf.gz")
	file, err := os.Create("testout.csv")
	if err != nil {
		fmt.Println(err)
		os.Exit(5)
	}
	writer := csv.NewWriter(file)
	defer writer.Flush()
	defer f.Close()
	defer gz.Close()
	//defer file.Close()
	//fmt.Println(sampleNames)
	//line, _ := reader.ReadString('\n')
	//line, _,_ := reader.ReadLine()
	//fmt.Println(reader.ReadString('\n'))
	//fmt.Println(reader.ReadString('\n'))
	for {
		line, err := reader.ReadString('\n')
		if err == io.EOF {
			break
		} else if err != nil {
			fmt.Println(err)
			os.Exit(1)
		}

		row := strings.Split(line, "\t")
		variant := Variant{
			Chr:    row[0],
			Pos:    row[1],
			Ref:    row[3],
			Alt:    row[4],
			Var_ID: "hi",
		}

		for i, _ := range row[9:] {
			//fmt.Println(i)
			call := strings.Split(row[i+9], ":")
			//fmt.Println(call)
			if call[0] == "0/1" {
				gt := strings.Split(call[2], ",")
				variant.Calls = append(variant.Calls, Call{
					PT_ID:  sampleNames[i],
					RNC:    call[7],
					GT:     call[0],
					GQ:     call[4],
					AD_Ref: gt[0],
					AD_Alt: gt[1],
				})
				//fmt.Println(variant.Calls)
			}
		}
		for i := range variant.Calls {
			writer.Write([]string{variant.Calls[i].PT_ID,
				variant.Calls[i].GT,
				variant.Calls[i].GQ,
				variant.Calls[i].AD_Ref,
				variant.Calls[i].AD_Alt,
				variant.Calls[i].RNC})
		}
	}

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
