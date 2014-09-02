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
import std.parallelism;
import std.process;
import std.c.string;


struct FastQRead {
	/*
		From bamread
	*/
	this(string read_name,                          // info for developers:
         string sequence,                           // these 3 fields are needed
         string quality_scores) {
        enforce(read_name.length < 256, "Too long read name, length must be <= 255");
        
        this._l_read_name = cast(ubyte)(read_name.length + 1); // tailing '\0'
        this._l_seq       = cast(int)(sequence.length);

        // set read_name
//        auto _offset = _read_name_offset;
//        _chunk[_offset .. _offset + read_name.length] = cast(ubyte[])read_name;
//        _chunk[_offset + read_name.length] = cast(ubyte)'\0';


        this.sequence = sequence;
	}
private:

    // Offsets of various arrays in bytes.
    // Currently, are computed each time, so if speed will be an issue,
    // they can be made fields instead of properties.
    @property size_t _read_name_offset() const nothrow pure { 
        return 8 * int.sizeof; 
    }

         
}

FastQRead parseAlignmentLine(string line, SamHeader header, OutBuffer buffer=null) {
	
	
}


FastQRead _parseFastQRecord(Tuple!(char[], FastQReader, OutBuffer) t) {
    auto r = parseAlignmentLine(cast(string)t[0], t[1]._header, t[2]);
    FastQRead result;
    if (t[1]._seqprocmode) {
        result = r;
    } else {
        auto storage = uninitializedArray!(ubyte[])(r.raw_data.length);
        storage[] = r.raw_data[];
        result.raw_data = storage;
    }
    result.associateWithReader(t[1]);
    return result;
}

void printUsage() {
    stderr.writeln("Usage: yamsvp-bedregion [options] <input.bam | input.sam>");
    stderr.writeln();
    stderr.writeln("Options: -T, --threads=NTHREADS");
    stderr.writeln("                    maximum number of threads to use");
    stderr.writeln("         -f, --fastq=sample.fastq");
    stderr.writeln("                    FastQ file to QC");
}


private {
    extern(C) size_t lseek(int, size_t, int);
    bool isSeekable(ref File file) {
        return lseek(file.fileno(), 0, 0) != ~0;
    }
}
class FastQReader {

    private {
        File openFastQFile(string filename) {
            if (filename.length < 4) {
                throw new Exception("invalid name for FastQ file: " ~ filename);
            } else if (filename[$ - 4 .. $] == ".bam") {
                throw new Exception("SAM reader can't read BAM file " ~ filename);
            } else {
                return File(filename);
            }
        }

    }

    ///
    this(string filename) {
        _file = openFastQFile(filename);
        _filename = filename;
        _seekable = _file.isSeekable();
        _initializeStream();
    }

    ///
    
    auto reads() @property {
    	_LineRange lines = _lines;
        if (_seekable) {
            if (_filename !is null) {
                auto file = openFastQFile(_filename);
                lines = ByLineFast(file);
            } else {
                _file.seek(0);
                lines = ByLineFast(_file);
            }
            auto dummy = lines.front;
            for (int i = 0; i < _lines_to_skip; i++)
                lines.popFront();
        }

        auto b = new OutBuffer(262144);
        return lines.zip(repeat(this), repeat(b)).map!_parseFastQRecord();
    }
    
    /// Filename
    string filename() @property const {
        return _filename;
    }

private:
    File _file;
    bool _seekable;
    string _filename;
    _LineRange _lines;
    ulong _lines_to_skip;

    bool _seqprocmode;

    void _initializeStream() {
        auto header = Appender!(char[])(); 

        _lines = ByLineFast(_file);

        while (!_lines.empty) {
            auto line = _lines.front;
            if (line.length > 0 && line[0] == '@') {
                header.put(line);
                header.put('\n');
                _lines_to_skip += 1;
                _lines.popFront();
            } else {
                break;
            }
        }
    }
}

int main(string[] args) {
	auto n_threads = 4;
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
    
    auto fastq = new FastQReader(args[1]);
    
    
}
