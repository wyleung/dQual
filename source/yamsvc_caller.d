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
    
    /*
    	read region definitions,
    	sort the regions by chromosome/contig, start
    		- filter by given DPdis
    		- filter out any clen < readlength, we are not capable calling these
    	chunk regions by: (n regions / n_treads) * 50?
    		- start analysis of each region in a new task
    		- tasks reports to RegionCollection class (shared ipc)
    	match regions
    	annotate regions (DEL, INS, INV, ITX, CTC etc.)
    	write vcf file
    */
    
	return 0;
}

