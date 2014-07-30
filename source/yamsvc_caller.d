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
import std.csv;
import std.file;

import bio.bam.referenceinfo;

import bio.bam.reader;
import bio.bam.pileup;
import bio.bam.read;
import bio.bam.writer;

import std.parallelism;

import yamsvc.utils;

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
    
    string bamfile = args[1];
    auto bam = new BamReader(bamfile, task_pool);
    
    // open a bed file
    auto bedrecords = readText(bed_file);
    auto records = csvReader!BedRecord(bedrecords, null , '\t');
    foreach(record; records){
        // select the discordant reads from this region
        
        CallRegion cr = new CallRegion();
        cr.bed = record;
        
        auto reads = bam[ record.chromosome ][cast(uint) record.start .. cast(uint) record.end];
        foreach( read; reads ) {
            if ( is_concordant(read) ){
                /* We don't want to analyse concordant reads, only discordant are used for SV calling */
                continue;
            }
            
            if ( !has_minimum_quality(read, 20)) {
                continue;
            }
            
            writefln("%s\t%d\t%s\t%d\tbin:%s\ttlen:%d\tflag:%d\tMQ:%d\tBaseCovered:%d", 
                read.ref_id, 
                read.position, 
                read.mate_ref_id, 
                read.mate_position, 
                read.bin, 
                read.template_length, 
                read.flag,
                read.mapping_quality, 
                read.basesCovered() );


            cr.adjust_pos( read.position, read.position+read.basesCovered() );
            
            cr.add_matechromosome( read.mate_ref_id );
            cr.add_mate_pos( read.mate_position );
            cr.addInsertsize( read.template_length );
        }
        writefln("%s\t%d\t%d\t%s\t%d",
            record.chromosome,
            record.start,
            record.end,
            record.infofields,
            cr.avgInsertsize
            );
    }
    
    
//    auto dst = new BamWriter(args[2], 9);     // maximal compression
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

