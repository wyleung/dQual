module yamsvc.utils;

import std.algorithm;
import std.stdio;
import std.parallelism;
import std.math;
import std.stream;

import bio.bam.read;
import bio.bam.writer;

import sambamba.utils.common.bed;
import sambamba.utils.common.filtering;




bool is_concordant( BamRead read ) {
	switch ( read.flag() ) {
		case 99:
		case 147:
		case 83:
		case 163:
			return true;
		default:
			return false;
	}
	return false;
}

bool has_minimum_quality( BamRead read, int minQuality ) {
	if(read.mapping_quality == 0) return false;
	return read.mapping_quality >= minQuality;
}

bool has_abnormal_isize( BamRead read, uint isize, uint stdev ) {
	if( ( abs(read.template_length) > ( isize + ( 2 * stdev ) ) ) || 
		( abs(read.template_length) < ( isize - ( 2 * stdev ) ) ) ) {
		return true;
	}
	return false;
}

bool isDiscordant( BamRead read, uint estimated_insertsize, uint estimated_insertsize_stdev ) {
	return (
	(
	 has_minimum_quality(read, 20) &&
	 read.is_paired &&
	 !(read.mate_is_unmapped) && 
	 (!is_concordant(read) || has_abnormal_isize(read, estimated_insertsize, estimated_insertsize_stdev)) )
	 );
}


// from sambamba/view.d (r4b7d32ba0de49f26d4ba0de8ba8b1e43c75bc50f)

void runProcessor(SB, R, F, P)(SB bam, R reads, F filter, P processor) {
    if (processor.is_serial)
        bam.assumeSequentialProcessing();
    if (cast(NullFilter) filter)
        processor.process(reads, bam);
    else
        processor.process(reads.filtered(filter), bam);
}

final class BamSerializer {
	/*
		This version is not Windows Compatible
	*/
    private string _f;
    private int _level;
    private TaskPool _task_pool;
    private enum BUFSIZE = 4096;//1_048_576;


    this(string f, int compression_level, TaskPool pool) {
        _f = f;
        _level = compression_level;
        _task_pool = pool;
    }

    enum is_serial = true;

    void process(R, SB)(R reads, SB bam) 
    {
        Stream output_stream = new BufferedFile(_f, FileMode.OutNew, 
                                                BUFSIZE);
        auto writer = new BamWriter(output_stream, _level, _task_pool);
        scope(exit) writer.finish();

        writer.writeSamHeader(bam.header);
        writer.writeReferenceSequenceInfo(bam.reference_sequences);
        foreach (read; reads)
            writer.writeRecord(read);
    }
}

import bio.bam.region;

BamRegion[] parseBed(Reader)(string bed_filename, Reader bam, bool non_overlapping=true) {
    auto index = sambamba.utils.common.bed.readIntervals(bed_filename, non_overlapping);
    BamRegion[] regions;
    foreach (reference, intervals; index) {
        if (!bam.hasReference(reference))
            continue;
        auto id = bam[reference].id;
        foreach (interval; intervals)
            regions ~= BamRegion(cast(uint)id,
                                 cast(uint)interval.beg, cast(uint)interval.end);
    }
    std.algorithm.sort(regions);
    return regions;
}


// flag filtering extended from sambamba


// Filtering integer fields
// based on the filters defined by sambamba 
final class AlignmentFlagFilter(string op) : Filter {
    private long _value;
    this(long value) {
        _value = value;
    }
    bool accepts(ref BamRead a) const {
        mixin("return a.flag " ~ op ~ "_value;");
    }
}

/// Filtering integer fields
final class TemplateSizeFilter(string op) : Filter {
    private long _value;
    this(long value) {
        _value = value;
    }
    bool accepts(ref BamRead a) const {
    	mixin("return std.math.abs(a.template_length) " ~ op ~ "_value;");
    }
}

struct Histogram {
public:
	void add( T )( T value ) {
		this.vals[ cast(uint) std.math.abs(value) ] += 1;
		this.total_count += 1;
	}
	
	@property uint[uint] raw() {
		return this.vals;
	}
	
    @property long size() {
        return this.total_count;
    }
    
    @property long width() {
        return this.vals.length;
    }
    
    @property double mean() {
        auto sum = 0;
        foreach(k,v; this.vals){
            sum += k*v;
        }
        return cast(double)(sum / this.size);
    }
    
    @property uint q25() {
        return this._q25;
    }
    
    @property uint q50() {
        return this._q50;
    }
    
    @property uint q75() {
        return this._q75;
    }
    
    @property double sd() {
        return this._sd;
    }
    
    @property uint min() {
        auto keys = this.vals.keys.sort;
        return keys[0];
    }
    @property uint max() {
        auto keys = this.vals.keys.sort;
        return keys[ keys.length - 1 ];
    }
    
    @property void flat_histogram() {
        auto keys = this.vals.keys.sort;
        foreach(k; keys){
        	writefln( "%s: %d", k, this.vals[k] );
    	}
	}
    
    uint[] _recompute() {
        auto q25_pos = this.size / 4;
        auto q50_pos = q25_pos*2;
        auto q75_pos = q25_pos*3;

        auto keys = this.vals.keys.sort;
        auto n_seen = 0;
        foreach( k; keys ) {
            foreach( val; 1 .. this.vals[k] ){
                ++n_seen;
//                writeln(n_seen);
//                writeln(q25_pos);
//                writeln(q50_pos);
//                writeln(q75_pos);
//                writeln(this._q75);

                if( n_seen == q25_pos ) {
                    this._q25 = k;
                }
                if( n_seen == q50_pos ) {
                    this._q50 = k;
                }
                if( n_seen == q75_pos ) {
                    this._q75 = k;
                }
            }
        }            
        this._sd = (this._q75 - this._q25) / 2 / 0.6744898;
        return [this._q25,this._q50,this._q75];
    }

	@property uint median() {
		auto mid = cast(uint) (this.size / 2);
		auto keys = this.vals.keys.sort;
		auto n_seen = 0;
		foreach( k; keys ) {
		    foreach( val; 1 .. this.vals[k] ){
		        mid--;
                if( mid == 0 ) {
                    this._median = k;
                    this._q50 = k;
                    return k;
                }
		    }
		}
		return 0;
	}
	
private:
	uint[ uint ] vals;
	uint total_count = 0;
    uint _q25 = 0;
    uint _q50 = 0;
    uint _q75 = 0;
    uint _mean = 0;
    uint _median = 0;
    double _sd = 0;
    
}
