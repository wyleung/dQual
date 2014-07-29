#!/usr/bin/env rdmd

/*
 @author: Wai Yi Leung
 @copyright: 2013-2014 Wai Yi Leung
 @contact: w.y.leung@e-sensei.nl
 @license: MIT Licence
*/

/*
    The YAMSVC caller in the D-language. Ported from the Python version.  
*/


import std.getopt;
import std.stdio;
import std.parallelism;
import std.algorithm : count, max, reduce;
import std.math : abs;

import bio.bam.referenceinfo;

import bio.bam.reader;
import bio.bam.pileup;
import bio.bam.read;
import bio.bam.writer;

import std.parallelism;

import yamsvc.utils: is_concordant, has_minimum_quality;

void printUsage() {
    stderr.writeln("Usage: yamsvp-caller [options] <input.bam | input.sam>");
    stderr.writeln();
    stderr.writeln("Options: -T, --threads=NTHREADS");
    stderr.writeln("                    maximum number of threads to use");
    stderr.writeln("         -b, --bedfile=BEDFILE");
    stderr.writeln("                    input bedfile created with yamsvc-region");
    stderr.writeln("         -v, --verbose");
    stderr.writeln("                    Turn on verbose mode");
}

class RegionDefiner {
    /**
    	Define a region from a alignment (bam file)
    
      Example:
      -------------------------------------------
      import std.parallelism, bio.bam.reader;
      void main() {
        auto pool = new TaskPool(4); // use 4 threads
        scope (exit) pool.finish();  // don't forget!
        auto region = new RegionDefiner("file.bam", pool);
        ...
      }
      -------------------------------------------
     */
	this(std.stream.Stream stream, 
         std.parallelism.TaskPool task_pool = std.parallelism.taskPool) {
        _task_pool = task_pool;
    }

    /// ditto
    this(string filename, std.parallelism.TaskPool task_pool) {
        _filename = filename;
//        _source_stream = getNativeEndianSourceStream();
//        this(_source_stream, task_pool);
    }

    /// ditto
    this(string filename) {
        this(filename, std.parallelism.taskPool);
    }

    /** Filename, if the object was created via file name constructor,
        $(D null) otherwise.
     */
    string filename() @property const {
        return _filename;
    }


    TaskPool _task_pool;
    size_t _buffer_size = 4096; // buffer size to be used for I/O

private:
    
    string _filename;                       // filename (if available)


}


void defineRegion(BamReader bam, ReferenceSequenceInfo reference_sequence) {
	auto reads = bam[ cast(string) reference_sequence.name() ][1 .. reference_sequence.length];
	auto pileup = makePileup(reads, true);
	foreach (column; pileup) {
		writefln("%s %s", column.coverage, reference_sequence.name());
		}
	}


void defineRegion(string bamfile, ReferenceSequenceInfo reference_sequence) {
	auto bam = new BamReader(bamfile);
	defineRegion(bam, reference_sequence);
	}



int main(string[] args) {
	int n_threads = 4;
	bool verbose;
	string bed_file;

    getopt(
        args,
        std.getopt.config.caseSensitive,
        "threads|T",  &n_threads,    // numeric
        "bedfile|b", &bed_file,
        "verbose|v", &verbose,   // flag
        );

    if (args.length < 2) {
        printUsage();
        return 0;
    }
    
    auto task_pool = new TaskPool(n_threads);
    scope(exit) task_pool.finish();
    
    auto bam = new BamReader(args[1], task_pool);
    string bamfile = args[1];
//    auto dst = new BamWriter(args[2], 9); 	// maximal compression
//    scope (exit) dst.finish();              // close the stream at exit
//    dst.writeSamHeader(bam.header);         // copy header and reference sequence info
//    dst.writeReferenceSequenceInfo(bam.reference_sequences);
    
    foreach (refseq; bam.reference_sequences) {
    	writefln("Starting thread for %s", refseq);
    	auto t = task!defineRegion(bamfile, refseq);
    	task_pool.put(t);
//    	t.executeInNewThread();
//    	defineRegion(bam, refseq);
    }
	return 0;
}

