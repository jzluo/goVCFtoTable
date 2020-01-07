package main

import (
	"bufio"
	"bytes"
	"compress/gzip"
	//"encoding/csv"
	"flag"
	"fmt"
	"io"
	"log"
	"os"
	"runtime"
	"runtime/pprof"
	//"strings"
	"time"
)

type Variant struct {
	Chr, Pos         []byte
	Ref, Alt, Var_ID []byte
	Calls            []Call
}

type Call struct {
	PT_ID, RNC             []byte
	GT, GQ, AD_Ref, AD_Alt []byte
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

//func ReadVCFgz(filename string) (*bufio.Reader, *os.File, *gzip.Reader, []string) {
func ReadVCFgz(filename string) (*bufio.Reader, *os.File, *gzip.Reader, [][]byte) {
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

	var sampleNames [][]byte

	for {
		//line, err := reader.ReadString('\n')
		line, err := reader.ReadBytes('\n')
		if err != nil {
			fmt.Println("Error parsing VCF file")
			fmt.Println(err)
			os.Exit(3)
		}

		//if strings.HasPrefix(line, "#CHROM") {
		if bytes.HasPrefix(line, []byte("#CHROM")) {
			sampleNames = ParseSampleNames(line)
			break
		}
	}

	return reader, f, gz, sampleNames
}

//func ParseSampleNames(line string) []string {
//var sampleNames []string
//splitLine := strings.Split(line, "\t")
//sampleNames = splitLine[9:]
//for i := range sampleNames {
//	sampleNames[i] = strings.Split(sampleNames[i], "_")[1]
//}
//return sampleNames
//}
func ParseSampleNames(line []byte) [][]byte {
	var sampleNames [][]byte
	splitLine := bytes.Split(line, []byte("\t"))
	sampleNames = splitLine[9:]
	for i := range sampleNames {
		sampleNames[i] = bytes.Split(sampleNames[i], []byte("_"))[1]
	}
	return sampleNames
}

func GetFields(call []byte) ([]byte, []byte, []byte, []byte) {
	//for {
	//	if bytes.IndexByte(call[4:], delim) == -1 {
	//
	//	}
	//	bytes.IndexByte(call[start])
	//}

	// get AD
	//fmt.Println("call", string(call))

	call = call[4:]
	end := bytes.IndexByte(call, ':')
	//end = bytes.IndexByte(call[4+start:], ',') + 4
	call = call[end+1:]
	end = bytes.IndexByte(call, ',')
	if end == -1 {
		return []byte("h"), []byte("h"), []byte("h"), []byte("h")
	}
	ref := call[:end]
	//fmt.Println("ref", string(ref))

	call = call[end+1:]
	end = bytes.IndexByte(call, ':')
	alt := call[:end]
	//fmt.Println("alt", string(alt))

	call = call[end+1:]
	//start = bytes.IndexByte(call[end:], ':') // SBPV end
	end = bytes.IndexByte(call, ':')
	call = call[end+1:] // skip SBPV
	end = bytes.IndexByte(call, ':')
	gq := call[:end]
	//fmt.Println("gq", string(gq))

	call = call[end+1:]
	end = bytes.IndexByte(call, ':')
	call = call[end+1:] // skip PL
	end = bytes.IndexByte(call, ':')
	call = call[end+1:]
	end = bytes.IndexByte(call, ':')
	call = call[end+1:] // skip FT
	//start = bytes.IndexByte(call[end:], ':') // PL end
	//end = bytes.IndexByte(call, ':')
	//start = bytes.IndexByte(call[end:], ':') // FT end
	rnc := call
	//fmt.Println("rnc", string(rnc))

	//t := bytes.Split(call, []byte(":"))
	//if len(t[2]) == 1 {
	//	return t[4], t[4], t[4], t[4]
	//}
	//ad := bytes.Split(t[2], []byte(","))
	////fmt.Println(string(call))
	//ref := ad[0]
	//alt := ad[1]
	//gq := t[4]
	//rnc := t[7]
	return ref, alt, gq, rnc
}

func SplitRow(row []byte) [][]byte {

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
	//writer := csv.NewWriter(file)
	writer := bufio.NewWriter(file)
	defer writer.Flush()
	defer f.Close()
	defer gz.Close()
	//defer file.Close()
	//fmt.Println(sampleNames)
	//line, _ := reader.ReadString('\n')
	//line, _,_ := reader.ReadLine()
	//fmt.Println(reader.ReadString('\n'))
	//fmt.Println(reader.ReadString('\n'))
	het := []byte("0/1")
	for {
		//line, err := reader.ReadString('\n')
		line, err := reader.ReadBytes('\n')
		if err == io.EOF {
			break
		} else if err != nil {
			fmt.Println(err)
			os.Exit(1)
		}

		//row := strings.Split(line, "\t")
		row := bytes.Split(line, []byte("\t"))
		variant := Variant{
			Chr:    row[0],
			Pos:    row[1],
			Ref:    row[3],
			Alt:    row[4],
			Var_ID: []byte("hi"),
		}

		for i, call := range row[9:] {
			//fmt.Println(i)
			//call := strings.Split(row[i+9], ":")
			//call := bytes.Split(row[i+9], []byte(":"))

			//fmt.Println(call)
			//if call[0] == "0/1" {
			//if bytes.Equal(call[0:3], []byte("0/1")) {
			if call[2] == 49 {
				ref, alt, gq, rnc := GetFields(call)
				//gt := strings.Split(call[2], ",")
				//gt := bytes.Split(call[2], []byte(","))
				variant.Calls = append(variant.Calls, Call{
					PT_ID:  sampleNames[i],
					RNC:    rnc,
					GT:     het,
					GQ:     gq,
					AD_Ref: ref,
					AD_Alt: alt,
				})
				//fmt.Println(variant.Calls)
			}
		}
		for i := range variant.Calls {
			//writer.Write([]string{variant.Calls[i].PT_ID,
			writer.Write(bytes.Join([][]byte{
				variant.Calls[i].GT,
				variant.Calls[i].GQ,
				variant.Calls[i].AD_Ref,
				variant.Calls[i].AD_Alt,
				variant.Calls[i].RNC,
			}, []byte(",")))
			writer.Write([]byte{'\n'})
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
