module yamsvc.datatypes;

import std.algorithm;
import std.math;
import std.conv;
import std.format;
import std.conv;
import std.range;
import std.stdio;
import std.string;
import std.format;

import bio.bam.read;

struct BedRecord {
	string chromosome;
	ulong start;
	ulong end;
	string infofields;
}

struct MatePair_t {
    BamRead * read1;
    BamRead * read2;
    
    this( BamRead * r ) {
        if ( (*r).is_first_of_pair ) {
            read1 = r;
        } else {
            read2 = r;
        }
    }
    
    bool overlaps( MatePair_t other ) {
        if ( this.overlaps_clusterA( other.read1 ) && this.overlaps_clusterB( other.read2 ) ) {
            return true;
        }
        return false;
    }

    bool overlaps_clusterA( BamRead * other ) {
        if ( this.read1.position <= other.position &&
            other.position <= (this.read1.position + this.read1.basesCovered) ) {
            return true;
        } else if (
            other.position <= (this.read1.position + this.read1.basesCovered)
            && 
            (this.read1.position + this.read1.basesCovered) <= (other.position + other.basesCovered)
            ) {
            return true;
        }
        return false;
    }
    bool overlaps_clusterB( BamRead * other ) {
        if ( this.read2.position <= other.position &&
            other.position <= (this.read2.position + this.read2.basesCovered) ) {
            return true;
        } else if (
            other.position <= (this.read1.position + this.read2.basesCovered)
            && 
            (this.read2.position + this.read1.basesCovered) <= (other.position + other.basesCovered)
            ) {
            return true;
        }
        return false;
    }

    bool opEquals()(auto ref const MatePair_t o ) const {
        // mate on overlap on both clusters
        // mate on overlap on insertsize with max 2 sd
        
    }
}

class DiscordantRegion {
	ulong start = 0;
	ulong end = 0;
	short chromosome=-1;
	ushort[] coverages_discordant;
	ushort[] coverages;
	bool working;
	ubyte orientation;
	
    uint[] mate_chromosome;
    ulong[] mate_positions;
    ushort[] orientations;
    uint[] insertsizes;

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

	ulong cov_disc;
	ulong cov_conc;
	
	this(ubyte _orientation){
		this.orientation = _orientation;
	}
	
	void addCoverage( ushort cov ) {
		this.coverages ~= cov;
		this.cov_conc += cov;
	}
	
	void addCoverageDiscordant( ushort cov ) {
		this.coverages_discordant ~= cov;
		this.cov_disc += cov;
	}
	
	ushort maxCoverage(){
		if(this.coverages.length){
			return reduce!(max)(this.coverages);
		} else {
			return 0;
		}
	}
	
	ushort maxCoverageDiscordant(){
		if(this.coverages_discordant.length){
			return reduce!(max)(this.coverages_discordant);
		} else {
			return 0;
		}
	}
	
	ushort avgCoverage(){
		/* FIXME: cast from ulong to int is not so performing? 
			We are only interested in coverages upto ushort (65,535) large?
		*/
		auto sum = reduce!((a,b) => a + b)(0, this.coverages);
		if( sum == 0 ) {
			return 0;
		}
		return cast(ushort) abs(sum / this.coverages.length);
	}
	
	@property ushort avgDisCoveragePerBase(){
		/* FIXME: cast from ulong to int is not so performing? 
			We are only interested in coverages upto ushort (65,535) large?
		*/
		auto sum = reduce!((a,b) => a + b)(0, this.coverages_discordant);
		if( sum == 0 ) {
			return 0;
		}
		return cast(ushort) abs(sum / this.size);
	}
	
	void resetCoverage() {
		/*
			Empty the coverages for a new region run
		*/
		this.coverages = [];
		this.coverages_discordant = [];
		this.cov_disc = 0;
		this.cov_conc = 0;
		this.mate_chromosome = [];
		this.mate_positions = [];
	}
	
	@property ulong size() {
		return this.end - this.start;
	}
	
	void startReporting( ushort ref_id, ulong start, ulong end ) {
		this.working = true;
		this.chromosome = ref_id;
		this.start = start;
		this.end = end;
		this.resetCoverage();
	}
	void startReporting( ) {
		this.working = true;
	}

    void toString(scope void delegate(const(char)[]) sink,
                  FormatSpec!char fmt)
    {
        put(sink, "[");
        formatValue(sink, this.start, fmt);
        put(sink, " : ");
        formatValue(sink, this.end, fmt);
        put(sink, ", orientation:");
        formatValue(sink, this.orientation, fmt);
        put(sink, ", ");
        put(sink, "]");
    }
}

unittest {
	auto dr = new DiscordantRegion(1);
	assert( dr.orientation == 1 );
	dr.addCoverage(1);
	assert( dr.cov_conc == 1 );
	dr.addCoverage(1);
	assert( dr.cov_conc == 2 );
	}


class CallRegion_old {
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



/**
    Lightweigth version of the $(bio.bam.read.BamRead),
    This version can be stored in memory with the essential information
*/

import std.bitmanip;

static struct ReadInfo {
    int ref_id; // max 4294967295 (4)
    uint pos; // max 4294967295 (4)
    ubyte basesCovered; // max 255 (1)

    int mate_ref_id; // max 4294967295 (4)
    uint mate_pos; // max 4294967295 (4)

    mixin(bitfields!(
        bool, "reverse",        1,
        bool, "mate_reverse",   1,
        bool, "first_of_pair",  1,
        uint, "",               5)); // .sizeof = ubyte

    ushort flag;
    string name;

//	~this() {
//		stderr.writeln("Del struct ReadInfo");
//	}

    @property uint start() {
        return this.pos;
    }

    @property uint end() {
        return this.pos+this.basesCovered;
    }

    @property bool transchromosomal() {
    	return this.ref_id != this.mate_ref_id;
	}

    this( BamRead r ) {
        this.name = r.name;
        this.ref_id = r.ref_id;
        this.pos = r.position;
        this.basesCovered = cast(ubyte) r.basesCovered;
        
        this.mate_ref_id = r.mate_ref_id;
        this.mate_pos = r.mate_position;
        
        this.reverse = r.is_reverse_strand;
        this.mate_reverse = r.mate_is_reverse_strand;
        this.first_of_pair = r.is_first_of_pair;
        
        this.flag = r.flag;
    }
//    
//    this( int ref_id, ubyte pos, ubyte basesCovered, ushort flag, string name ) {
//        this.ref_id = ref_id;
//        this.pos = pos;
//        this.basesCovered = basesCovered;
//        this.flag = flag;
//        this.name = name;
//    }
}

