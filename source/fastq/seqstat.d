#!/usr/bin/env rdmd

/*
 @author: Wai Yi Leung
 @copyright: 2013-2014 Wai Yi Leung
 @contact: w.y.leung@e-sensei.nl
 @license: MIT Licence
*/
import std.array;
import std.getopt;
import std.stdio;
import std.parallelism;
import std.algorithm : count, max, reduce;

import bio.core.utils.outbuffer;
import bio.core.utils.bylinefast;
alias ByLineFast _LineRange;


import std.string;
import std.range;
import std.algorithm;
import std.typecons;
import std.process;
import std.c.string;
import std.json;

class SeqStat {

    uint[string] bases;
    uint[string] reads;
    uint[string] length;
    
    this() {
        length["len_min"] = uint.max;
        length["len_max"] = uint.min;
        reads["N"] = 0;
        reads["total"] = 0;
        bases["N"] = 0;
        bases["total"] = 0;
    }
    
    void add( FastQRead read ) {
        // record n bases
        this.bases["N"] += read.n_bases;
        this.bases["total"] += read.length;
        
        if( read.n_bases > 0 ){
            reads["N"] += 1;
        }
        this.reads["total"] += 1;
        
        // update length statistics
        if ( read.length < this.length["len_min"] ) {
            this.length["len_min"] = read.length;
        }
        if ( read.length > this.length["len_max"] ) {
            this.length["len_max"] = read.length;
        }
        
    }
    
    void report() {
//    	JSONValue j = JSONValue();
//    	j["stats"] = cast(JSONValue) "";
//    	writeln(j);
        writefln( "%s \t %s \t %s\n", this.bases, this.reads, this.length );
        }
    
}


class FastQRead {
    /*
        From bamread
    */
    immutable string readname;
    string sequence;
    immutable string holder;
    string quality_scores;
    
    this(string readname,                          // info for developers:
         string sequence,                           // these 3 fields are needed
         string holder,
         string quality_scores) {
//        enforce(read_name.length < 256, "Too long read name, length must be <= 255");

        this.readname = readname;
        this.sequence = sequence;
        this.holder = holder;
        this.quality_scores = quality_scores;
    }
    
    @property ushort avgqual() {
        auto sum = reduce!((a,b) => a + b)(0, this.quality_scores);
        if( sum == 0 ) {
            return 0;
        }
        return cast(ushort) (sum / this.quality_scores.length);

    }
    
    @property uint n_bases() {
        auto _n = 0;
        foreach( base; this.sequence ) {
                if( base == 'N' ) {
                    _n += 1;
                }
            }
        return _n;
        }
    
    @property uint[] quals() {
        uint[] ret;
        
        foreach(qual; this.quality_scores){
            ret ~= cast(uint) qual;
        }
        return ret;
    }
    
    @property uint length() {
        return cast(uint) this.sequence.length;
        }
    @property string record() {
        return format("%s\n%s\n%s\n%s\n", this.readname, this.sequence, this.holder, this.quality_scores);
        }
}

void printUsage() {
    stderr.writeln("Usage: yamsvp-bedregion [options] <input.fastq>");
    stderr.writeln();
    stderr.writeln("Options: -T, --threads=NTHREADS");
    stderr.writeln("                    maximum number of threads to use");
    stderr.writeln("         -f, --fastq=sample.fastq");
    stderr.writeln("                    FastQ file to QC");
}

void doStats( SeqStat stats, ref FastQRead read ) {
    stats.add(read);
    }

int main(string[] args) {
    
    auto n_threads = 1; // actually any number higher than 1 is now dangerous, skipping tasks etc. DEBUG!
    string s_fastq;
    
    getopt(
        args,
        std.getopt.config.caseSensitive,
        "threads|T",  &n_threads,    // numeric
        "fastq|F",  &s_fastq,    // numeric
        );

    if (args.length < 2) {
        printUsage();
        return 0;
    }
    
    auto task_pool = new TaskPool(n_threads);
    scope(exit) task_pool.finish();
    
    File file = File(args[1], "r");
        
    SeqStat stats = new SeqStat();

    string record;
    char[] buf;
    
    while (file.readln(buf)) {
        string name = chomp(cast(string)buf);
        string sequence = chomp(file.readln());
        string holder = chomp(file.readln());
        string quality = chomp(file.readln());
        
        auto r = new FastQRead(name, sequence, holder, quality);
        doStats( stats, r );
//        auto t = task!doStats( stats, r );
//        task_pool.put(t);
//      write(r.record);
//      writefln("Read length: %d", r.length);
//      writefln("Base qualities: %s", r.quals);
//      writefln("Avg. quality: %d", r.avgqual);
//      writefln("N bases: %d", r.n_bases);
        
    }
    stats.report();
    
    return 0;
}
