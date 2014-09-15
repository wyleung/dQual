module bam.SortedBamReader;

import bio.bam.read;
import bio.bam.reader;
import bio.bam.readrange;
import bio.bam.reference;
import bio.bam.randomaccessmanager;
import bio.bam.region;

import std.stdio;
import std.range;
import std.exception;

struct AlnPair_t {
	/*
	The convention read numbering:
	We follow the ordering from the sequencer, so the first read in pair is read1 and the second read2
	*/
    BamRead read1;
    BamRead read2;
    
    this (BamRead r1, BamRead r2) {
    	this.read1 = r1;
    	this.read2 = r2;
    }
    
    @property string name() const {
    	return read1.name;
	}
    
    @property bool is_duplicate() const {
    	return read1.is_duplicate || read2.is_duplicate;
    }
    
    @property bool is_secondary_alignment() const {
    	return read1.is_secondary_alignment || read2.is_secondary_alignment;
    }

    @property bool is_unmapped() const {
    	return read1.is_unmapped || read2.is_unmapped;
    }
    
    @property bool is_interchromosomal() const {
    	return read1.ref_id != read2.ref_id;
    }
    
    @property uint insertsize() const {
    	// the insert size definition the following pair configurations:
    	// no insertsize on interchromosomal matings
    	if( this.is_interchromosomal ) {
    		return 0;
		}
    	
    	// F1 F2: F2.start - F1.start
    	if( !read1.is_reverse_strand && !read2.is_reverse_strand ) {
    		if ( read1.position > read2.position ) {
    			return (read1.position + read1.basesCovered) - read2.position;
			}
    		return read2.position - read1.position;
		}
    	
    	// R1 R2: R1.end - R2.end
    	if( read1.is_reverse_strand && read2.is_reverse_strand ) {
    		if ( read2.position > read1.position ) {
    			return (read2.position + read2.basesCovered) - read1.position;
			}
    		return (read1.position+read1.basesCovered) - (read2.position+read2.basesCovered);
		}
    	
    	// F1 R2: R2.end - F1.start
    	if( !read1.is_reverse_strand && read2.is_reverse_strand ) {
    		if ( read1.position > read2.position ) {
    			return (read1.position + read1.basesCovered) - read2.position;
			}
    		return (read2.position + read2.basesCovered) - read1.position;
		}
    	// F2 R1: R1.end - F2.start
    	if( read1.is_reverse_strand && !read2.is_reverse_strand ) {
    		if ( read2.position > read1.position ) {
    			return (read2.position + read2.basesCovered) - read1.position;
			}
    		return (read1.position + read1.basesCovered) - read2.position;
		}
    	return 0;
    }
    
}

struct ReadPairRange {
	/* support on non sliced bam file (e.g. no chromosome selected) 
	* we work directly on the reader.reads
	*/
	
    ulong reads_processed = 0;

	this(BamReader reader=null) {
        _reader = reader;
        _reads = reader.reads;
        readNext();
    }

    @property bool empty() const {
        return this._paired_reads.length == 0;
    }

    @property ref AlnPair_t front() {
    	BamRead r1 = this._paired_reads[0];
    	BamRead r2 = this._paired_reads[1];
    	
		this.currentpair = AlnPair_t(r1, r2);
    	return this.currentpair;
    }
    
//    @property ref BamRead front() {
//    	return this._paired_reads[0];
//    }
    
    void popFront() {
		this._paired_reads = this._paired_reads[2 .. $];
		/*
			We can read ahead for caching purpose, upto 2000 pairs in the cache
		*/
    	if ( !this._reads.empty && this._paired_reads.length < 4000 ) {
    		this.readNext();
		}
    }
    
    void readNext() {
        this._reads.popFront();
        this._current_record = this._reads.front.dup; // have a dup, so that no reference is stored to the last read?

        bool found = false;
        int i = 0;
        foreach(int j, BamRead cread; this._cached_reads) {
            if( cread.name == this._current_record.name) {
                found = true;
                this._cached_reads = std.algorithm.remove( this._cached_reads, j );
                this._paired_reads ~= cread;
                this._paired_reads ~= this._current_record;
            }
        } /* end foreach */
        if ( !found ) {
            // if we have reached this part, then the current_read (this._read) was not found in the _cached_reads
            // append the read to the _cached_reads;
            this._cached_reads ~= this._current_record.dup;
            // recursively call readNext till we have a pair?
            this.readNext();
        }
    }
    
    /* end range functions */
    
private:
    BamReader _reader;
    BamRead _current_record;
    BamReadRange!withoutOffsets _reads;
    BamRead[] _cached_reads;
    BamRead[] _paired_reads;
    AlnPair_t currentpair;
} 

class SortedBamReader {
	/// currently not taking chromosome slices because of incompatible datatype (BamReadRange!withoutOffsets)
private:
    BamReader _bamreader;
public:
    // constructor
    this( string filename ) {
        BamReader reader = new BamReader(filename, std.parallelism.taskPool);
        this( reader );
    }

    /* ditto */
    this( BamReader reader ) {
        this._bamreader = reader;
    }

    /* bamfile helper functions */
    bool has_index() {
    	return this._bamreader.has_index();
	}

	auto opIndex(string ref_name) {
        enforce(this._bamreader.hasReference(ref_name), "Reference with name " ~ ref_name ~ " does not exist");
        return 0;
//        auto ref_id = _reference_sequence_dict[ref_name];
//        return reference(ref_id);
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
}

