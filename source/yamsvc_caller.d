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

//import core.memory;
import std.algorithm;
import std.concurrency;
import std.csv;
import std.datetime;
import std.getopt;
import std.math : abs;
import std.parallelism;
import std.path;
import std.stdio;
import std.string;

//import std.file;
//import std.stream;

import bio.bam.pileup;
import bio.bam.read;
import bio.bam.reader;
import bio.bam.referenceinfo;
import bio.bam.writer;

import sambamba.utils.common.filtering;

import yamsvc.utils;
import yamsvc.datatypes;

void printUsage() {
    stderr.writeln("Usage: yamsvc-caller [options] <input.bam | input.sam>");
    stderr.writeln();
    stderr.writeln("Options: -T, --threads=NTHREADS");
    stderr.writeln("                    maximum number of threads to use");
    stderr.writeln("         -o, --bam_output=BAMOUTPUT");
    stderr.writeln("                    Write the resulting bam to ...");
    stderr.writeln("         -b, --bedfile=BEDFILE");
    stderr.writeln("                    input bedfile created with yamsvc-region");
    stderr.writeln("         -l, --compression_level");
    stderr.writeln("                    compression level 1-9 (low-high) of bam file [9]");
    stderr.writeln("         -i, --insertsize");
    stderr.writeln("                    insertsize of the input bam [450]");
    stderr.writeln("         -s, --insertsd");
    stderr.writeln("                    stdev in the insertsize [15]");
    stderr.writeln("         -v, --verbose");
    stderr.writeln("                    Turn on verbose mode");
}

uint[] extractReads( string bamfile, BedRecord bed, 
uint estimated_insertsize, uint estimated_insertsize_stdev, 
TaskPool taskpool, BamWriter bamout ) {
    auto tp = new TaskPool(4);
    scope(exit) tp.finish();
    
auto bam = new BamReader(bamfile, tp);

    auto input_buf_size = 64_000_000;
    bam.setBufferSize(input_buf_size);
//bam.assumeSequentialProcessing();

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
    ubyte compression_level = 9;
    bool verbose;
    string bam_outputfilename;
    string bed_filename;

    getopt(
        args,
        std.getopt.config.caseSensitive,
        "threads|T",  &n_threads,    // numeric
        "bedfile|b", &bed_filename,
        "bamoutput|o", &bam_outputfilename,
        "compression_level|l", &compression_level, // numeric
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
//
//    auto writer_pool = new TaskPool(n_threads);
//    scope(exit) writer_pool.finish();
    
    string bamfile = args[1];
    BamReader bam = new BamReader(bamfile, reader_pool);
    // the buf-size is in bytes, 64 Kb = 1024 * 1024 * 64
auto input_buf_size = 32_000_000;
bam.setBufferSize(input_buf_size);
    
    /*
    ulong output_buf_size = 32_000_000;
    Stream stream = new BufferedFile(bam_output, FileMode.OutNew, 
                                         output_buf_size);
    BamWriter bamout = new BamWriter(stream, 9, writer_pool, output_buf_size); // maximal compression
    scope (exit) bamout.finish();              // close the stream at exit
    bamout.writeSamHeader(bam.header);         // copy header and reference sequence info
    bamout.writeReferenceSequenceInfo(bam.reference_sequences);
    */
//
//    File bam_outputfile; 
//    if (bam_outputfilename is null)
//        bam_outputfile = stdout;
//    else
//        bam_outputfile = File(bam_outputfilename, "w+");
//    

    Filter read_filter = new NullFilter();
    read_filter = new AndFilter(read_filter, new ValidAlignmentFilter());
    
    Filter only_mapped = new AndFilter(
new NotFilter(new FlagFilter!"is_unmapped"()), 
new NotFilter(new FlagFilter!"mate_is_unmapped"())
);
    
    Filter quality = new IntegerFieldFilter!">="("mapping_quality", 20);
    
    Filter flags = new NotFilter(
    new OrFilter( new AlignmentFlagFilter!"=="(99), 
    new OrFilter( new AlignmentFlagFilter!"=="(147), 
    new OrFilter( new AlignmentFlagFilter!"=="(83), 
    new AlignmentFlagFilter!"=="(163) ) ) )
);
    
    Filter template_size = new OrFilter( 
new TemplateSizeFilter!">="( (estimated_insertsize + ( 2 * estimated_insertsize_stdev )) ),
new TemplateSizeFilter!"<="( (estimated_insertsize - ( 2 * estimated_insertsize_stdev )) )
);
    
    
    read_filter = new AndFilter(read_filter, only_mapped);
    read_filter = new AndFilter(read_filter, quality);
    read_filter = new AndFilter(read_filter, new OrFilter( flags, template_size )); 
    
    /* This query is to extract the reads for downstream processing with the 
    Python based caller. The selection  and export is not fully needed in the D implementation.
    
    */

    int processAlignments(P)(P processor) {
    
        void runProcessor(SB, R, F)(SB bam, R reads, F filter) {
            if (processor.is_serial)
                bam.assumeSequentialProcessing();
    writeln("Before processing ");
            if (cast(NullFilter) filter)
                processor.process(reads, bam);
            else
                processor.process(reads.filtered(filter), bam);
    writeln("After processing ");
        }
        
    // adapted from sambamba
    auto regions = parseBed(bed_filename, bam);
    auto reads = bam.getReadsOverlapping(regions);
    writeln("After reads");
    
runProcessor(bam, reads, read_filter);
    writeln("Before return");
return 0;
    }

writeln("here");
int stat =processAlignments(new BamSerializer(bam_outputfilename, cast(int) compression_level, reader_pool));
writeln("before exit");
return stat;

       
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
}

