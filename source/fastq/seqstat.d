#!/usr/bin/env rdmd

/*
 @author: Wai Yi Leung
 @copyright: 2013-2014 Wai Yi Leung
 @contact: w.y.leung@e-sensei.nl
 @license: MIT Licence
*/
import core.atomic;
import core.thread;
import core.time;

import std.algorithm : map, count, max, reduce, filter, equal;
import std.array;
import std.c.string;
import std.conv;
import std.digest.sha;
import std.format;
import std.json;
import std.getopt;
import std.math;
import std.process;
import std.stdio;
import std.string;
import std.typecons;
import std.range;
import std.parallelism;

import bio.core.utils.outbuffer;
import bio.core.utils.bylinefast;
alias ByLineFast _LineRange;

import yamsvc.gzip;
import yamsvc.utils: Histogram;

class SeqStat {

    private ulong[string] bases;
    private ulong[string] reads;
    private ulong[string] length;
    private uint[8] base_quality; // representing classes 1, 10 - 70
    private uint[8] read_quality; // representing classes 1, 10 - 70
    private Histogram histo_read = Histogram();
    private Histogram histo_base = Histogram();
    private uint phred_correction = 33;
    private string path;
    private string sha1sum;

    this() {
        length["len_min"] = ulong.max;
        length["len_max"] = ulong.min;
        reads["N"] = 0;
        reads["total"] = 0;
        bases["N"] = 0;
        bases["total"] = 0;
    }
    ~this() {
        writeln("Killed seqstat!");
        }
    
    public void setPath( string s ) {
    	this.path = s;
	}

    public void setSHA1sum( string s ) {
    	this.sha1sum = toLower(s);
	}
    
    public string phred_encoding() {
    	uint l_qual = this.histo_base.min;
    	uint h_qual = this.histo_base.max;

    	// first determine whether this is sanger
    	if( h_qual > 74 ) {
    		// could be one of the solexa, illumina 1.3+, illumina 1.5+
    		this.phred_correction = 64; // correct the offset
    		return "solexa";
    	} else if ( l_qual < 59 ) {
    		// we have a 'sanger' type of sequence
    		this.phred_correction = 33;
    		return "sanger";
    	} 
    	
    	return "Unknown";
    }
    
    public void add( FastQRead read ) {
        // record n bases
        foreach( i; 0 .. read.n_bases ) {
            this.bases["N"]++;
        }
        foreach( i; 0 .. read.length ) {
            this.bases["total"]++;
        }
        // keep track of the base quality add them to a histogram
        // used to fill this.base_quality after the last read (where also determine the phred encoding)
        foreach( qual; read.quals ) {
        	uint _q = cast(uint) abs(qual);
        	this.histo_base.add( _q );
        }
        this.histo_read.add( read.avgqual );

        reads["total"]++;
        if( read.n_bases > 0 ){
            reads["N"]++;
        }
        
        // update length statistics
        if ( read.length < this.length["len_min"] ) {
            this.length["len_min"] = read.length;
        }
        if ( read.length > this.length["len_max"] ) {
            this.length["len_max"] = read.length;
        }
        
    }
    
    public void report() {
    	// get encoding, and calculate offset:
    	string encoding = this.phred_encoding;
    	auto base_keys = this.histo_base.raw.keys.sort;
    	
        foreach(k; base_keys){
        	uint phred_qual = k - this.phred_correction;
//        	writefln( "%d: %d", phred_qual, this.histo_base.raw[k]);

	        foreach( i; 0 .. this.base_quality.length ) {
	            auto qval = i * 10u ? i * 10u : 1u;  // to fix the 0 index,resulting into qval=0
	            if( phred_qual >= qval ) {
	            	this.base_quality[ i ] += this.histo_base.raw[ k ];
	            }
            }
    	}
        
    	auto read_keys = this.histo_read.raw.keys.sort;
    	
        foreach(k; read_keys){
        	uint phred_qual = k - this.phred_correction;
//        	writefln( "%d: %d", phred_qual, this.histo_read.raw[k]);

	        foreach( i; 0 .. this.read_quality.length ) {
	            auto qval = i * 10u ? i * 10u : 1u;  // to fix the 0 index,resulting into qval=0
	            if( phred_qual >= qval ) {
	            	this.read_quality[ i ] += this.histo_read.raw[ k ];
	            }
            }
    	}
//    	writefln("base raw counts: %s", this.histo_base.raw);
//    	writefln("base histogram width: %d", this.histo_base.width);
//    	writefln("base min, max, phred-encoding: %d, %d, %s, -%d", 
//    		this.histo_base.min, this.histo_base.max, this.phred_encoding, phred_correction);
//    	
//    	writefln("bases counted: %d", this.histo_base.size);
//        writefln( "%s \t %s \t %s \t %s\n", this.bases, this.base_quality, this.reads, this.length );

        auto base_quals = appender!string();
        base_quals.put("{\n");
        
        foreach( i; 0 .. this.base_quality.length ) {
            auto qval = i * 10u ? i * 10u : 1u;
        	base_quals.put( "\"" );
        	base_quals.put( to!string(qval) );
        	base_quals.put( "\"" );
        	base_quals.put( ": " );
        	base_quals.put( to!string(this.base_quality[i]) );
        	if( i == this.base_quality.length-1 ) {
        		base_quals.put( "\n" );
        	}
        	else{
        		base_quals.put( ",\n" );
        	}
        }
        base_quals.put("}");
        
        auto read_quals = appender!string();
        read_quals.put("{\n");
        
        foreach( i; 0 .. this.read_quality.length ) {
            auto qval = i * 10u ? i * 10u : 1u;
        	read_quals.put( "\"" );
        	read_quals.put( to!string(qval) );
        	read_quals.put( "\"" );
        	read_quals.put( ": " );
        	read_quals.put( to!string(this.read_quality[i]) );
        	if( i == this.read_quality.length-1 ) {
        		read_quals.put( "\n" );
        	}
        	else{
        		read_quals.put( ",\n" );
        	}
        }
        read_quals.put("}");
        
        auto writer = appender!string();
        
        formattedWrite(writer, "{
	\"files\": {
		\"fastq\": {
			\"checksum_sha1\": \"%s\",
			\"path\": \"%s\"
		}
	},
	\"stats\": {
		\"qual_encoding\": \"%s\",
		\"bases\": {
			\"num_n\": %d,
			\"num_qual_gte\": %s,
			\"num_total\": %d
		},
		\"reads\": {
			\"len_max\": %d,
			\"len_min\": %d,
			\"num_mean_qual_qte\": %s,
			\"num_total\": %d,
			\"num_with_n\": %d
		}
	}
}", 
		this.sha1sum,
		this.path,
		this.phred_encoding,
        this.bases["N"],
		base_quals.data,
        this.bases["total"],
        this.length["len_max"],
        this.length["len_min"],
        read_quals.data,
        this.reads["total"],
        this.reads["N"]
        );
        writeln(writer.data);

        //		JSONValue myjson = parseJSON( writer.data );
//		writeln(toJSON(&myjson, true));   
    }
}

final class FastQRead {
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
        foreach(qual; this.quality_scores) {
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
    stderr.writeln("Usage: fastq-seqstat [options] <input.fastq>");
    stderr.writeln();
    stderr.writeln("Options: -T, --threads=NTHREADS");
    stderr.writeln("                    maximum number of threads to use");
    stderr.writeln();
}

void printBytes(T)(ref T variable)
{
    const ubyte * begin = cast(ubyte*)&variable;    // (1)

    writefln("type   : %s", T.stringof);
    writefln("value  : %s", variable);
    writefln("address: %s", begin);                 // (2)
    writef  ("bytes  : ");

    writefln("%(%02x %)", begin[0 .. T.sizeof]);    // (3)

    writeln();
}
class Lock
{}
void printMemory(T)(T * location, size_t length)
{
    const ubyte * begin = cast(ubyte*)location;

    foreach (address; begin .. begin + length) {
        char c = (isPrintable(*address) ? *address : '.');

        writefln("%s:  %02x  %s", address, *address, c);
    }
}

void countReads(FastQRead read, shared(uint) * counter, SeqStat stats, shared(Lock) lock) {
    atomicOp!"+="(*counter, 1);    // atomic update
    stats.add(read);
//    synchronized(lock){
//        atomicOp!"+="(*counter, 1);    // atomic update
//    }
}

void doStats( SeqStat stats, FastQRead read, shared(int) * counter, shared(Lock) lock ) {
    
    synchronized (stats) {
        stats.add(read);
    }
    
    printBytes(*counter);
    atomicOp!"+="(*counter, 1);    // atomic update
    synchronized (lock) {
        ++(*counter);
    }
    printBytes(*counter);
}

void makeReport(SeqStat *stats) {
    stats.report();
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
    
    string input_file = args[1];
    std.file.File file;
    auto len = input_file.length;
    if(input_file[len - 3 .. len] == ".gz") {
        auto pipe = pipeShell("gunzip -c " ~ input_file);
        file = pipe.stdout;
//    	file = GzipByLine(input_file);
	} else {
    	file = File(input_file, "r");
	}
    

    SeqStat stats = new SeqStat();

    string record;
    char[] buf;
    
    shared(uint) reads_processed = uint.min;
    shared(uint) * ptr = &reads_processed;
    shared(Lock) lock = new shared(Lock)();
    
    SHA1 sha;
    sha.start();
    string line;
    while ( (line = file.readln()) !is null ) {
        string name = line;
        string sequence = file.readln();
        string holder = file.readln();
        string quality = file.readln();
        
        auto read = new FastQRead(
            chomp(name), 
            chomp(sequence), 
            chomp(holder),
            chomp(quality));
        countReads( read, ptr, stats, lock );

        sha.put(cast(immutable(ubyte)[]) name);
        sha.put(cast(immutable(ubyte)[]) sequence);
        sha.put(cast(immutable(ubyte)[]) holder);
        sha.put(cast(immutable(ubyte)[]) quality);;
    }

    task_pool.finish();
    stats.setPath(args[1]);
    stats.setSHA1sum(toHexString(sha.finish()));
    stats.report();
    return 0;
}
