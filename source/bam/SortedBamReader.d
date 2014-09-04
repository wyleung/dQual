module bam.SortedBamReader;

import bio.bam.read;
import bio.bam.reader;
import bio.bam.readrange;
import bio.bam.reference;
import bio.bam.randomaccessmanager;
import bio.bam.region;

import std.stdio;

struct AlnPair_t {
    BamRead read1;
    BamRead read2;
//    ReferenceSequenceInfo refinfo;
//    string name;
}

struct ReadPairRange {
    ulong reads_processed = 0;

	this(BamReader reader=null) {
        _reader = reader;
        _reads = reader.reads;
//        writefln("Here");
        // initialize by filling the first pair info.
        readNext();
//        writefln("After");
    }

    @property bool empty() const {
//        writefln("empty");
        return this._paired_reads.length == 0;
    }
    
    @property ref BamRead front() {
//        writefln("front");
//        return this._current_record;
    	return this._paired_reads[0];
    	
    }
    
    void popFront() {
//        writefln("popFront");
		this._paired_reads = this._paired_reads[1 .. $];
		/*
			We can read ahead for caching purpose, upto 4 pairs in the cache
		*/
    	if ( !this._reads.empty && this._paired_reads.length < 8 ) {
    		this.readNext();
		}
    }
    
    void readNext() {
    	reads_processed++;
//    	writeln(reads_processed);
        this._reads.popFront();
        this._current_record = this._reads.front.dup; // have a dup, so that no reference is stored to the last read?

//        writefln( "Record: %s", this._current_record.toSam() );

        bool found = false;
        int i = 0;
        foreach(int j, BamRead cread; this._cached_reads) {
            if( cread.name == this._current_record.name) {
//            	writefln("Found: %s", cread.name);
                found = true;
                this._cached_reads = std.algorithm.remove( this._cached_reads, j );
                this._paired_reads ~= cread;
                this._paired_reads ~= this._current_record;
            }
        } /* end foreach */
        if ( !found ) {
            // if we have reached this part, then the current_read (this._read) was not found in the _cached_reads
            // append the read to the _cached_reads;
//            writefln( "Not found: %s, adding to cache", this._current_record.name );
            this._cached_reads ~= this._current_record.dup;
//            writefln( "%d in cache ", this._cached_reads.length );
            // recursively call readNext till we have a pair?
            this.readNext();
        }
    }
    
    /* end range functions */
    
private:
    BamReader _reader;
    BamRead _current_record;
    BamReadRange!withoutOffsets _reads;
    bool _empty = false;
    BamRead[] _cached_reads;
    BamRead[] _paired_reads;
} 

class SortedBamReader {
	/// currently not taking chromosome slices because of incompatible datatype (BamReadRange!withoutOffsets)
private:
    BamReader _bamreader;
    BamRead _read;
    BamReadRange!withoutOffsets _reads;
//    ReferenceSequence _reads;
    BamRead[] _cached_reads;
public:
    // constructor
    this( string filename ) {
        BamReader reader = new BamReader(filename, std.parallelism.taskPool);
        this( reader );
    }
    
    /* ditto */
//    this( string filename, string chromosome ) {
//        BamReader reader = new BamReader(filename, std.parallelism.taskPool);
//        ReferenceSequenceInfo ref_id;
//        foreach( _ref_id; reader.reference_sequences ) {
//            if( _ref_id.name == chromosome ) {
//                ref_id = _ref_id;
//                break;
//            }
//        }
//        this( reader, ref_id ); 
//    }
    
    /* ditto */
    this( BamReader reader ) {
        this._bamreader = reader;
        this._reads = reader.reads();
    }
    
    /* ditto */
//    this( BamReader reader, ReferenceSequenceInfo ref_id ) {
//        this._bamreader = reader;
//        this._reads = this._bamreader[ ref_id.name ][0 .. ref_id.length];
//    }

    /* bamfile helper functions */
    bool has_index() {
    	return this._bamreader.has_index();
	}
    
    void createIndex() {
    	this._bamreader.createIndex();
	}
    /* end bamfile helper functions */
    

    auto reads() @property {
//    	return this._bamreader.reads;
    	return ReadPairRange( this._bamreader );
	}

    const(bio.bam.referenceinfo.ReferenceSequenceInfo)[] reference_sequences() @property const {
        return this._bamreader.reference_sequences;
    }
    
    AlnPair_t processNextRead() {
        while(1) {
            bool found_read = false;
            foreach(int j, BamRead cread; this._cached_reads) {
                if( cread.name == this._read.name) {
                    found_read = true;
                    this._cached_reads = std.algorithm.remove( this._cached_reads, j );
                    AlnPair_t newPair = AlnPair_t( cread, this._read );
                    return newPair;
                }
            } /* end foreach */
            
            // if we have reached this part, then the current_read (this._read) was not found in the _cached_reads
            // append the read to the _cached_reads;
            
            this._cached_reads ~= this._read.dup;
            
//            this.advance();
        }
    }
}

