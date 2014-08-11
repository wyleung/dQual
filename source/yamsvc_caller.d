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
import std.datetime;
import std.stdio;
import std.parallelism;
import std.algorithm;
import std.math : abs;
import std.csv;
import std.file;
import std.concurrency;
import std.path;
import std.file;
import std.stream;
import std.stdio;

import bio.bam.referenceinfo;

import bio.bam.reader;
import bio.bam.pileup;
import bio.bam.read;
import bio.bam.writer;

import yamsvc.utils;

void printUsage() {
    stderr.writeln("Usage: yamsvp-caller [options] <input.bam | input.sam>");
    stderr.writeln();
    stderr.writeln("Options: -T, --threads=NTHREADS");
    stderr.writeln("                    maximum number of threads to use");
    stderr.writeln("         -b, --bedfile=BEDFILE");
    stderr.writeln("                    input bedfile created with yamsvc-region");
    stderr.writeln("         -i, --insertsize");
    stderr.writeln("                    insertsize of the input bam [450]");
    stderr.writeln("         -s, --insertsd");
    stderr.writeln("                    stdev in the insertsize [15]");
    stderr.writeln("         -v, --verbose");
    stderr.writeln("                    Turn on verbose mode");
}

class CallRegion {
    BedRecord bed;
    
    uint[] mate_chromosome;
    ulong[] mate_positions;
    uint[] insertsizes;
    uint clusterID;
    ushort[] orientations;
    uint start;
    uint end;
    
    uint[][] adj_pos;
    
    /* not necessary needed, only when this is re-evaluated for the 2nd time */
    void adjust_pos ( uint s_pos, uint e_pos ) {
        this.adj_pos ~= [s_pos, e_pos]; 
    }
    
    void addOrientation( short orientation ) {
        this.orientations ~= orientation;
    }
    
    void add_matechromosome( uint matechr ) {
        this.mate_chromosome ~= matechr;
    }
    void add_mate_pos( ulong matepos ) {
        this.mate_positions ~= matepos;
    }

    void addInsertsize( uint insertsize ) {
        this.insertsizes ~= abs(insertsize);
    }

    @property uint avgInsertsize(){
        auto sum = reduce!((a,b) => a + b)(0, this.insertsizes);
        if( sum == 0 ) {
            return 0;
        }
        return cast(uint) abs(sum / this.insertsizes.length);
    }
}

uint[] extractReads( string bamfile, BedRecord bed, 
					uint estimated_insertsize, uint estimated_insertsize_stdev, 
					TaskPool taskpool, BamWriter bamout ) {
    auto tp = new TaskPool(4);
    scope(exit) tp.finish();
    
	auto bam = new BamReader(bamfile, tp);

	auto input_buf_size = 64_000_000;
	bam.setBufferSize(input_buf_size);
//	bam.assumeSequentialProcessing();

	auto reads = bam[ bed.chromosome ][cast(uint) bed.start .. cast(uint) bed.end];
    auto reads_processed = 0;
    auto reads_read = 0;

    BamRead[] _reads;

    foreach( read; reads ) {
    	reads_read++;
	    if ( !isDiscordant(read, estimated_insertsize, estimated_insertsize_stdev ) ){
	        /* We don't want to analyse concordant reads, only discordant are used for SV calling */
	        continue;
	    }
		reads_processed++;
		_reads ~= read;
    }
    foreach( r; _reads ){
	    bamout.writeRecord( r );
    }
    
    return [reads_processed, reads_read];
}

int main(string[] args) {
    int n_threads = 4;
    uint estimated_insertsize = 450;
    uint estimated_insertsize_stdev = 15;
    bool verbose;
    string bed_file;

    getopt(
        args,
        std.getopt.config.caseSensitive,
        "threads|T",  &n_threads,    // numeric
        "bedfile|b", &bed_file,
        "insertsize|i", &estimated_insertsize, // numeric
        "insertsd|s", &estimated_insertsize_stdev, // numeric
        "verbose|v", &verbose,   // flag
        );

    if (args.length < 2) {
        printUsage();
        return 0;
    }
    
    auto reader_pool = new TaskPool(n_threads);
    scope(exit) reader_pool.finish();

    auto writer_pool = new TaskPool(n_threads);
    scope(exit) writer_pool.finish();
    
    string bamfile = args[1];
    BamReader bam = new BamReader(bamfile, reader_pool);
    
    
    ulong output_buf_size = 32_000_000;
    Stream stream = new BufferedFile(args[2], FileMode.OutNew, 
                                         output_buf_size);
    BamWriter bamout = new BamWriter(stream, 9, writer_pool, output_buf_size); // maximal compression
    scope (exit) bamout.finish();              // close the stream at exit
    bamout.writeSamHeader(bam.header);         // copy header and reference sequence info
    bamout.writeReferenceSequenceInfo(bam.reference_sequences);

    // open a bed file
    auto bedrecords = readText(bed_file);
    auto records = csvReader!BedRecord(bedrecords, null , '\t');
    
    auto reads_processed = 0;
    auto reads_read = 0;
    auto records_processed = 0;
    
    
    stderr.writeln("Looping over bed records");
    StopWatch sw;
    sw.start();
    
    foreach(record; records){
        // select the discordant reads from this region
        if (record.end <= record.start) {
        	auto e = record.start;
        	record.start = record.end;
        	record.end = e+1;
        } 
        auto stats = extractReads( bamfile, 
					record, 
					estimated_insertsize, 
					estimated_insertsize_stdev, 
					reader_pool, 
					bamout );
        
        reads_processed += stats[0];
        reads_read += stats[1];
        /*
        auto t = task!extractReads( bamfile, 
									record, 
									estimated_insertsize, 
									estimated_insertsize_stdev, 
									reader_pool, 
									bamout );
        reader_pool.put(t);
        */
        CallRegion cr = new CallRegion();
        cr.bed = record;
    	records_processed++;
    	
    	if ( records_processed % 500 == 0 ) {
    		auto time_used = sw.peek().msecs;
    		writefln( "Processed %d bed records, exported %d reads in %d seconds, avg. %d read/sec, read %d reads/msec",
    			 records_processed, reads_processed, time_used/1000, reads_processed*1000/time_used, reads_read/time_used);
    	}
    }
    
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

