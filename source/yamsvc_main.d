#!/usr/bin/env rdmd

/*
 @author: Wai Yi Leung
 @copyright: 2013-2014 Wai Yi Leung
 @contact: w.y.leung@e-sensei.nl
 @license: MIT Licence
*/

/*
    Region definition: a genomic region with at least 2 bases
    We summarize the genome regions which have genomic DNA reads aligned, 
    but they are reported discordant (abnormal insert-size, mate on other chromosome)
*/


import core.atomic;
import core.memory;

import std.array;
import std.algorithm;
import std.datetime;
import std.file;
import std.getopt;
import std.math;
import std.parallelism;
import std.path;
import std.stdio;
import std.string;

import bio.bam.pileup;
import bio.bam.read;
import bio.bam.reader;
import bio.bam.region;
import bio.bam.writer;

import sambamba.utils.common.filtering;

import bam.SortedBamReader;
import yamsvc.datatypes;
import yamsvc.utils;

struct DiscordantCall {
    int ref1;
    uint s1; // cluster 1 start
    uint e1; // cluster 1 end

    int ref2;
    uint s2; // ditto
    uint e2; // ditto
    
    // translation table for orientation:
    // 0: FF
    // 1: FR
    // 2: RF
    // 3: RR
    ubyte orientation;
    ushort pairs_supporting;
    static const ushort SEARCH_WINDOW = 500;
    
    string repr() {
        auto ret = appender!string();
        ret.put( to!string(ref1) ); ret.put(" ");
        ret.put( to!string(this.s1) ); ret.put(" ");
        ret.put( to!string(this.e1) ); ret.put(" ");
        ret.put( to!string(ref2) ); ret.put(" ");
        ret.put( to!string(this.s2) ); ret.put(" ");
        ret.put( to!string(this.e2) ); ret.put(" ");
        ret.put("- orientation: ");
        ret.put( to!string(orientation) ); ret.put(" ");
        ret.put("- insertsize: ");
        ret.put( to!string(this.insertsize) ); ret.put(" ");
        ret.put("- DP: ");
        ret.put( to!string( this.pairs_supporting ) ); ret.put(" ");
        return cast(string) ret.data;
    }
    
    @property string toBED( const(ReferenceSequenceInfo)[] refseq ) const {
        return format("%s\t%d\t%d\tDP=%d;SVLEN=%d\n%s\t%d\t%d\tDP=%d;SVLEN=%d", 
            refseq[ref1].name, s1, e1, pairs_supporting, insertsize, 
            refseq[ref2].name, s2, e2, pairs_supporting, insertsize);
    }
    
    this( AlnPair_t pair ) {
        this.ref1 = pair.read1.ref_id;
        this.s1 = pair.read1.position;
        this.e1 = pair.read1.position + pair.read1.basesCovered;
        
        this.ref2 = pair.read2.ref_id;
        this.s2 = pair.read2.position;
        this.e2 = pair.read2.position + pair.read2.basesCovered;
        
        this.orientation = pair.orientation;
        this.pairs_supporting = 1;
    }
    
    @property uint scan_till_pos() {
        return this.e2 + this.SEARCH_WINDOW;
    }
    
    bool overlapping( BamRead other, uint reference_id, uint start, uint end ) const {
        if ( other.ref_id != reference_id ) return false;
        // compute the start and end of the 'other' first
        uint os = other.position;
        uint oe = other.position + other.basesCovered;
        
        /** /overlapping by start/
         *          sccccccccccccccceWWWWWWW <- cluster
         *                    srrrrrrrrrrrrrre <- read
         */
        bool overlapping_by_start = start <= os && os <= (end + this.SEARCH_WINDOW);
        
        /** /overlapping by end/
         *          sccccccccccccccceWWWWWWW <- cluster
         *     srrrrrrrrrrrrrrre
         ****/
        bool overlapping_by_end = oe >= start && oe <= end + this.SEARCH_WINDOW;
        bool enclosed = start <= os  && oe <= (end + this.SEARCH_WINDOW);
        return overlapping_by_start || overlapping_by_end || enclosed; 
    }
    
    
    
    bool overlaps( AlnPair_t pair ) const {
        if ( pair.orientation != this.orientation ) return false;
//        writeln(pair.repr());
        
        double size_ratio;
        if ( this.insertsize == 0 ) {
        	size_ratio = cast(double)(pair.insertsize / (this.insertsize+1));
    	} else {
    		size_ratio = cast(double)(pair.insertsize / this.insertsize);
		}
        
        
        return ( size_ratio <= 1.1 && size_ratio >= 0.9 ) && 
                overlapping( pair.read1, ref1, s1, e1 ) && 
                overlapping( pair.read2, ref2, s2, e2 );
    }

    void merge( AlnPair_t pair ) {
        if ( pair.read1.position < s1 ) this.s1 = pair.read1.position;
        if ( pair.read2.position < s2 ) this.s2 = pair.read2.position;

        uint newE1 = pair.read1.position + pair.read1.basesCovered;
        uint newE2 = pair.read2.position + pair.read2.basesCovered;
        if ( newE1 > e1 ) this.e1 = newE1;
        if ( newE2 > e2 ) this.e2 = newE2;
        this.pairs_supporting += 1;
    }

    bool opEquals()(auto const ref DiscordantCall other) const {
        return ( other.ref1 == this.ref1 ) &&
                ( other.s1 == this.s1 ) &&
                ( other.e1 == this.e1 ) &&
                ( other.ref2 == this.ref2 ) &&
                ( other.s2 == this.s2 ) &&
                ( other.e2 == this.e2 ) &&
                ( other.orientation == this.orientation );
    }

    int opCmp( ref const DiscordantCall other ) const {
        if ( this == other ) {
            return 0;
        } else if ( this.ref1 < other.ref1 ) {
            return -1;
        } else if ( this.ref2 < other.ref2 ) {
            return -1;
        } else if ( this.s1 < other.s1 ) {
            return -1;
//        } else if ( this.e1 < other.e1 ) {
//            return -1;
//        } else if ( this.s2 < other.s2 ) {
//            return -1;
//        } else if ( this.e2 < other.e2 ) {
//            return -1;
        } else  {
            return 1;
        }
    }
    
    @property bool is_interchromosomal() const {
        return ref1 != ref2;
    }
    
    @property int insertsize() const {
        // the insert size definition the following pair configurations:
        // no insertsize on interchromosomal matings
        if( this.is_interchromosomal ) {
            return 0;
        }
        
        switch( orientation ) {
            default:
                return 0;
            case 0: // FF
                return s2 - s1;
            case 1: // FR
                return e2 - s1;
            case 2: // RF
                return s2 - e1;
            case 3: // RR
                return e2 - e1;
        }
    }
}

unittest {
    // testing the DiscordantCall
    // first create 3 AlnPair_t to use
}

struct ReadStorage {
private:
	BamRead[] _reads;
public:
    @property bool empty() const {
        return this._reads.length == 0;
    }

    @property ref BamRead front() {
    	return this._reads[0];
    }
   
    void popFront() {
        this._reads = this._reads[1 .. $];
    }
    
    void opOpAssign(string op)(BamRead r)
	{
	    mixin("_reads" ~ op ~ "=r;");
	}
	
	void add( BamRead r ) {
		try {
			this._reads ~= r.dup;
		}
		catch (Exception e) {
			stderr.writefln("[Error] %s", e.msg);
		}
	}
	
}

ReadInfo fromBamRead( BamRead r ) {
    ReadInfo ri = ReadInfo( r );
    return ri;
}

void extractReadInfo( string bamfile, 
                string discordant_bam, // target bam
                string calls_file,  // target calls file
                string bed_file,  // target bed file
                BamRegion bamregion, 
                ushort cutoff_min, 
                ushort cutoff_max, 
                ushort window_width,
                uint estimated_insertsize, 
                uint estimated_insertsize_stdev,
                uint max_sd_degree) {
//    stderr.writefln("Entering extract read info %s", bamregion);
    BamReader bam;
    auto tp = new TaskPool(4);
    scope(exit) tp.finish();

    try {
        bam = new BamReader(bamfile, tp);
    } catch (Exception e) {
        writefln("[Error] %s", e.msg);
    }
    auto input_buf_size = 16_000_000;
    bam.setBufferSize(input_buf_size);
    auto reads = bam.getReadsOverlapping([bamregion]);

    BamRead[] contigreads_trans;
    BamRead[] contigreads_first;
    BamRead[] contigreads_second;
    BamRead read;
    AlnPair_t pair;
    DiscordantCall[] clusters;
    DiscordantCall[] clusters_for_reporting;
    
    // resource allocation in memory
    static const ushort CAPACITY_SIZE = 10240;
    contigreads_trans.reserve(CAPACITY_SIZE * 4);
    contigreads_first.reserve(CAPACITY_SIZE);
    contigreads_second.reserve(CAPACITY_SIZE);
    clusters.reserve(10240);
    clusters_for_reporting.reserve(10240);
    
    uint reads_seen = 0;
//    defaultPoolThreads(1);
    std.datetime.StopWatch sw;
    sw.start(); // total time spent on reading all the reads in the region
    auto _minireads = tp.map!fromBamRead(
            reads
        ).array();
    auto trans_reads = filter!( a => a.transchromosomal == true )(_minireads).array();
    auto minireads = filter!( a => a.transchromosomal == false )(_minireads).array();
    sw.stop();
    
    
    std.datetime.StopWatch swsorting;
    swsorting.start(); // total time spent on reading all the reads in the region

    auto readsFirst = filter!( a => a.first_of_pair == true )(minireads).array();
    sort!("a.pos < b.pos")(readsFirst);

    auto readsSecond = filter!( a => a.first_of_pair == false )(minireads).array();
    sort!("a.mate_pos < b.mate_pos")(readsSecond);

    stderr.writefln("Sorting:\t%s\treads in total %d [hns],\t%g [hns]/read",
        bamregion, 
        swsorting.peek().hnsecs, 
        cast(double) swsorting.peek().hnsecs / minireads.length);

    stderr.writefln("Mapped:\t%s\t:%d reads in total %d [hns],\t%g [hns]/read\nFirst:\t%d\tSecond\t%d",
        bamregion, 
        minireads.length, 
        sw.peek().hnsecs, 
        cast(double) sw.peek().hnsecs / minireads.length,
        readsFirst.length,
        readsSecond.length );
    
    ReadInfo[] report_storage;
    report_storage.reserve(10240*2);
    
    ReadInfo[] tmp_storage;
    tmp_storage.reserve(10240*2);
    ReadInfo last_read;
	ReadInfo r1;
	ReadInfo r2;
	AlnPair_t[] cluster_pairs;
	cluster_pairs.reserve(10240);
	uint cluster_count;

    std.datetime.StopWatch sw_matching;
    sw_matching.start(); // total time spent on reading all the reads in the region

    stderr.writefln("1: %d, 2: %d", readsFirst.length, readsSecond.length );

    while ( !readsFirst.empty && !readsSecond.empty ) {
    	r1 = readsFirst.front;
    	readsFirst.popFront();

    	if ( r1.pos != last_read.pos ) {
    		// dump the tmp_storage to 
    		report_storage ~= tmp_storage;
    		writefln("cur: %d\tlast: %d", r1.pos, last_read.pos);
    		foreach( r; tmp_storage) {
    		    writeln(r);
    		    }
    		tmp_storage.length = 0;
    		tmp_storage.reserve(1024);
    		writeln("Restarting tmp storage");

    		// fill the tmp_storage will all reads second in pair on the particular position of the r1.
    		// a one time operation
            // assume the readsSecond is sorted by mate_pos
            while ( readsSecond.front.mate_pos == r1.pos ) {
                tmp_storage ~= readsSecond.front;
                readsSecond.popFront();
            }
		}
        last_read = r1;

    	if ( tmp_storage.length == 0 ) {
    	    // no mating reads for this position, move on to next?
    	    report_storage ~= r1;
    	    continue;
	    }


    	auto hits = filter!( a => a.name == r1.name )( tmp_storage ).array();
    	if( hits.length ) {
    		r2 = hits.front;
    		writefln("%s\t%s", r1.name, r2.name);
//    		cluster_pairs ~= AlnPair_t(r1, r2);
    		cluster_count++;
    		tmp_storage = filter!( a => a.name != r1.name )( tmp_storage ).array();
		} else {
		    writefln("%d\t%d", tmp_storage.length, hits.length);
//		    tmp_storage ~= r1;
	    }
	}
    stderr.writefln("Matching:\t%s\t%d pairs in total %d [hns],\t%g [hns]/pair",
        bamregion, 
        cluster_count,
        sw_matching.peek().hnsecs, 
        cast(double) sw_matching.peek().hnsecs / cluster_count);
    
    

    return; // quick hack to not run the code below.




















    while (!reads.empty) {
        reads_seen++;
        reads.popFront();
        read = reads.front.dup; // duplicate only for low memory usage ?, dereferencing the read, optimize for the popFront
        bool found_read = false;

        if ( read.ref_id != read.mate_ref_id ) {
            if ( contigreads_trans.length == contigreads_trans.capacity ) {
                contigreads_trans.reserve( contigreads_trans.length + CAPACITY_SIZE );
            }
            contigreads_trans ~= read;
            continue;
        }

//	    std.datetime.StopWatch sw_matching;
//	    sw_matching.start(); // total time spent on reading all the reads in the region

        BamRead[] workset;
        if( read.is_first_of_pair ) {
            workset = contigreads_second;
        } else {
            workset = contigreads_first;
        }
        
        uint readposition = read.position;
        string readname = read.name;
        auto matching = filter!( a => 
        						a.mate_position == readposition &&
        						a.name == readname )(workset);

        foreach ( BamRead match; matching ) {
            found_read = true;
            pair = AlnPair_t(match, read);
            break;
        }
//	    sw_matching.stop();
//	    stderr.writefln("%s : %d reads in total %d [hns],  %g [hns]/read",
//	        bamregion, 
//	        workset.length, 
//	        sw_matching.peek().hnsecs, 
//	        cast(double) sw_matching.peek().hnsecs / workset.length );

        //        if (matching.length > 0) {
//            found_read = true;
//            pair = AlnPair_t(matching.front, read);
//        }
//        foreach(int j, BamRead cread; contigreads) {
//            if( cread.mate_position != position || cread.name != readname ) continue;
//            found_read = true;
//            pair = AlnPair_t(cread, read);
//            contigreads = std.algorithm.remove( contigreads, j );
//            break;
//        }

        if (!found_read) {
            if( read.is_first_of_pair ) {
                if ( contigreads_first.length == contigreads_first.capacity ) {
                	stderr.writefln("1: Capacity old: %d", contigreads_first.capacity);
                	contigreads_first.reserve( contigreads_first.length + CAPACITY_SIZE );
                	stderr.writefln("1: Capacity new: %d", contigreads_first.capacity);
            	}
                contigreads_first ~= read;
            } else {
                if (contigreads_second.length == contigreads_second.capacity ) {
                	stderr.writefln("2: Capacity old: %d", contigreads_second.capacity);
                	contigreads_second.reserve( contigreads_second.length + CAPACITY_SIZE );
                	stderr.writefln("2: Capacity new: %d", contigreads_second.capacity);
            	}
                contigreads_second ~= read;
            }
        } else { // the read was found
    	}
        
    	if ( reads_seen % CAPACITY_SIZE == 0 ) {
    		// cleanup the read storage every 10k reads processed.
    		// using the position of the current read to determine which ones to keep
//    		ulong now = sw.peek().hnsecs;
//    		ulong f_seq = contigreads_first.length;
//    		ulong s_seq = contigreads_second.length;
//    		stderr.writefln("%d\tReads: \t%d \t1: \t%d \t%d \t2: \t%d \t%d",
//    		    read.position,
//    			reads_seen,
//    			contigreads_first.length, contigreads_first.capacity,
//    			contigreads_second.length, contigreads_second.capacity
//    			);
    		
            contigreads_first = sort!("a.mate_position < b.mate_position")(
            						filter!( a => a.mate_position >= read.position )(contigreads_first).array()
        						).array();
        	contigreads_first.reserve( contigreads_first.length + CAPACITY_SIZE );
            
            
            contigreads_second = sort!("a.mate_position < b.mate_position")(
            						filter!( a => a.mate_position >= read.position )(contigreads_second).array()
        						).array();
            contigreads_second.reserve( contigreads_second.length + CAPACITY_SIZE );
            
            
//            ulong time_lapsed = (sw.peek().hnsecs - now);
//    		stderr.writefln("%d\tReads: \t%d \t1: \t%d \t%d \t2: \t%d \t%d\t\t%d\t%d\t%d",
//    		    read.position,
//    			reads_seen,
//    			contigreads_first.length, contigreads_first.capacity,
//    			contigreads_second.length, contigreads_second.capacity,
//    			time_lapsed,
//    			f_seq - contigreads_first.length,
//    			s_seq - contigreads_second.length
//    			);
        }
        
//        else {
//            bool report = false;
//            if( pair.is_unmapped  || pair.is_duplicate || pair.is_secondary_alignment ) continue;
//
//            // don't use translocational mates for insertsize histogram
//            if( pair.is_interchromosomal ) report=true;
//
//            uint isize_low = estimated_insertsize - (max_sd_degree * estimated_insertsize_stdev );
//            uint isize_high = estimated_insertsize + (max_sd_degree * estimated_insertsize_stdev );
//            if( report || pair.insertsize <= isize_low || pair.insertsize >= isize_high  ) {
//                bool found_pair = false;
//                int[] to_remove;
//                version(Debug) {
//                    writefln("Scanning cluster (%d) matching the pair %s", clusters.length, pair.name);
//                }
//                // remove anything that is not in the scanning window, only if the number of elements is exceeding n=5000
//                if ( clusters.length >= 5000 ) {
//                    clusters_for_reporting ~= clusters.filter!( a => a.scan_till_pos < pair.read1.position ).array();
//                    clusters = clusters.filter!( a => a.scan_till_pos >= pair.read1.position ).array();
//                    sort!("a < b")(clusters);
//                }
//                std.datetime.StopWatch sw_clustering;
//                sw_clustering.start(); // measure time spent on seeking matching pair
//
//                foreach( int callno, ref DiscordantCall call; clusters ) {
//                    if( !call.overlaps( pair ) ) continue;
//
//                    found_pair = true;
//                    call.merge( pair );
//                    break;
//                }
//                sw_clustering.stop();
//                writefln("Time iterating over clusters: %d [hns] and %s", sw_clustering.peek().hnsecs, found_pair);
//                
//                if ( !found_pair ) {
//                    clusters ~= DiscordantCall( pair );
//                } // if not found pair
//            } //if in insertsize range or reporting
//        } // found read, else ...
    } // end while

    sw.stop();
    stderr.writefln("%s : %d reads in total %d [hns],  %g [hns]/read",
        bamregion, 
        reads_seen, 
        sw.peek().hnsecs, 
        cast(double) sw.peek().hnsecs / reads_seen );
    stdout.writefln("%s : %d reads in total %d [hns],  %g [hns]/read",
        bamregion, 
        reads_seen, 
        sw.peek().hnsecs, 
        cast(double) sw.peek().hnsecs / reads_seen );

    clusters ~= clusters_for_reporting;

    auto f = File(calls_file, "w"); // open for writing
    auto b = File(bed_file, "w"); // open for writing
    
    sort!("a < b")(clusters);
    auto filtered_clusters = filter!(x => x.pairs_supporting >= 4)(clusters);
    foreach( DiscordantCall call; filtered_clusters) {
    		f.writefln("%s", call.repr());
    		b.writefln("%s", call.toBED( bam.reference_sequences ));
    }

    auto dst = new BamWriter(discordant_bam, 9); // maximal compression
    scope (exit) dst.finish();              // close the stream at exit
    dst.writeSamHeader(bam.header);         // copy header and reference sequence info
    dst.writeReferenceSequenceInfo(bam.reference_sequences);

//    writefln( "Unmatched reads: %d in %s, size in mem: %d", contigreads.length, reference_sequence.name(), contigreads.sizeof*contigreads.length );
    
    stderr.writefln("%d transchromosomal reads", contigreads_trans.length);
    stdout.writefln("%d transchromosomal reads", contigreads_trans.length);
    foreach( BamRead _r; contigreads_trans ) {
        dst.writeRecord( _r );
    }
    foreach( BamRead _r; contigreads_first ) {
        dst.writeRecord( _r );
    }
    foreach( BamRead _r; contigreads_second ) {
        dst.writeRecord( _r );
    }

}

void scanContig(string bamfile,
                string workdir,
                ReferenceSequenceInfo reference_sequence,
                uint refno,
                ushort cutoff_min, 
                ushort cutoff_max, 
                ushort window_width,
                uint binsize,
                uint estimated_insertsize, 
                uint estimated_insertsize_stdev,
                uint max_sd_degree,
                uint threads) {
    
//    stderr.writefln("Starting analysis of %s", reference_sequence.name);
//    stderr.writefln("Read bam %s for %s", bamfile, reference_sequence.name);

    // we create 2 files in the workdir:
    // new bam with only the discordant reads (for debugging)
    // CallRegion file with the plain calls without normalisation and filtering
    // just compose the filename over here and pass them to the analysis module.

    string cleaned_ref_name = replace(reference_sequence.name,"|","_");
    string wd = absolutePath(workdir ~ "/" ~ cleaned_ref_name ~ "/"); 
    // first make the workdir
    if ( !std.file.exists(wd) ) {
        try {
            mkdirRecurse( wd );
        }
        catch  (Exception e) {
//            fail silently, since the directory could be created by other scanContig threads
//            stderr.writefln("[Error] %s", e.msg);
        }
    }
    uint searchbins = reference_sequence.length / binsize;
//    auto tp = new TaskPool(1);
//    scope(exit) tp.finish();

//    stderr.writefln("[Debug] Chromosome %s len: %d, bins: %d, maxsize: %d", reference_sequence.name, reference_sequence.length, searchbins, (searchbins+1)*binsize);

    foreach( uint bin; 0 .. searchbins+1 ) {
        uint start = 0 + (bin*binsize);
        uint end = start + binsize;
        if (end > reference_sequence.length) {
            end = reference_sequence.length;
        }
        
//        stderr.writefln("%s %d %d %d", reference_sequence.name, reference_sequence.length, start, end);
        BamRegion bamregion = BamRegion( refno, start, end );
//        stderr.writefln("Bam region: %s", bamregion);
        string discordant_bam = absolutePath( wd ~ "/" ~ cleaned_ref_name ~ "." ~ to!string(bin) ~ ".bam" );
        string calls = absolutePath( wd ~ "/" ~ cleaned_ref_name ~ "." ~ to!string(bin) ~ ".calls" );
        string bedfile = absolutePath( wd ~ "/" ~ cleaned_ref_name ~ "." ~ to!string(bin) ~ ".bed" );

//        auto t = task!
        extractReadInfo(bamfile, 
                discordant_bam,
                calls,
                bedfile,
                bamregion, 
                cutoff_min, 
                cutoff_max, 
                window_width,
                estimated_insertsize, 
                estimated_insertsize_stdev,
                max_sd_degree);
//        tp.put(t);
    }

//    writefln("Writing new bam to %s", discordant_bam);
}

void printUsage() {
    int n_threads = std.parallelism.totalCPUs;
    
    stderr.writeln("Usage: yamsvp-bedregion [options] <input.bam | input.sam>");
    stderr.writeln();
    stderr.writeln("Options: -T, --threads=NTHREADS");
    stderr.writefln("                    maximum number of threads to use [%d]", n_threads);
    
    stderr.writeln("         -w, --window=WINDOWWIDTH");
    stderr.writeln("                    width of window (bp) determing overlap [200]");
    stderr.writeln("         -b, --binsize=BINSIZE");
    stderr.writeln("                    width of bin (bp) to scan [10_000_000]");

    stderr.writeln("         -l, --cutoff_min");
    stderr.writeln("                    low threshold for reads covering a breakpoint [2]");
    stderr.writeln("         -h, --cutoff_max");
    stderr.writeln("                    high threshold for reads covering a breakpoint [90]");

//    stderr.writeln("         -o, --bam_output=BAMOUTPUT");
//    stderr.writeln("                    Write the resulting bam to ...");
    stderr.writeln("         -l, --compression_level");
    stderr.writeln("                    compression level 1-9 (low-high) of bam file [9]");

    stderr.writeln("         -r, --reads_for_histogram");
    stderr.writeln("                    reads used to determine insert-size, a minimum of 10_000 is recomended [10_000]");
    stderr.writeln("         -d, --max_sd_degree");
    stderr.writeln("                    set degree of sd allowed to mark as concordant read in the insertsize [2]");

    stderr.writeln("         -v, --verbose");
    stderr.writeln("                    Turn on verbose mode");
}

int main(string[] args) {

    int n_threads = std.parallelism.totalCPUs;
    uint binsize = 10_000_000;
    ushort window_width = 200;

    ushort cutoff_max = 90;
    ushort cutoff_min = 2;

    uint estimated_insertsize = 450;
    uint estimated_insertsize_stdev = 15;
    uint max_sd_degree = 4;
    ubyte compression_level = 9;
    
    uint reads_for_histogram = 10_000;
    
    string bam_outputfilename;
    string bed_filename;

    bool verbose;
    string workdir = "./work/";

    getopt(
        args,
        std.getopt.config.caseSensitive,
        "threads|T",  &n_threads,    // numeric

        "window|w", &window_width,   // numeric
        "binsize|b", &binsize,   // numeric

        "cutoff_max|h", &cutoff_max, // numeric
        "cutoff_min|l", &cutoff_min, // numeric

        "bamoutput|o", &bam_outputfilename,
        "compression_level|l", &compression_level, // numeric

        "insertsize|i", &estimated_insertsize, // numeric
        "insertsd|s", &estimated_insertsize_stdev, // numeric
        "reads_for_histogram|r", &reads_for_histogram, // numeric
        "max_sd_degree|d", &max_sd_degree, // numeric

        "workdir|W", &workdir,       // flag
        "verbose|v", &verbose,       // flag
        );

    if (args.length < 2) {
        printUsage();
        return 0;
    }

//    std.parallelism.defaultPoolThreads = n_threads;
//    int working_threads = n_threads;
    int working_threads = n_threads;
//	if ( n_threads > 1 ) {
//		working_threads = cast(int)std.math.ceil( cast(float)(n_threads/2.0) );
//	}
    auto task_pool = new TaskPool( working_threads );
    scope(exit) task_pool.finish();
    
    string _bamfilename = args[1];
    auto bam = new BamReader(_bamfilename );
    auto src = new SortedBamReader( bam );
    
    // Fix index, we need this for random access
    if( !src.has_index ) {
        stderr.writeln("No bai index found, creating one now....");
        src.createIndex();
    }
    
    
    Histogram h = Histogram();
    auto _reads = src.reads;
    auto limit = reads_for_histogram;
    int counter = 0;
    foreach( AlnPair_t pair; _reads) {
    	++counter;
    	
    	if( pair.is_duplicate ) continue;
        if( pair.is_secondary_alignment) continue;
        if( pair.is_unmapped ) continue;
    	
    	// don't use translocational mates for insertsize histogram
        if( pair.is_interchromosomal ) continue;

        limit--;
        h.add( std.math.abs( pair.insertsize ) );
        
//        h.add( std.math.abs(read.template_length) );
        if (limit == 0){
            break;
        }
    }
    
    auto stat = h._recompute();
    
    stderr.writefln("Histogram range: %d .. %d", h.min , h.max);
    stderr.writefln("Histogram width: %d", h.width);
    stderr.writefln("Histogram records: %d", h.size);
    stderr.writefln("Histogram median: %s", h.median);
    stderr.writefln("Histogram q25: %s", h.q25);
    stderr.writefln("Histogram q50: %s", h.q50);
    stderr.writefln("Histogram q75: %s", h.q75);
    stderr.writefln("Histogram mean: %s", h.mean);
    stderr.writefln("Histogram sd: %s", h.sd);
    
    ReferenceSequenceInfo[] ref_seqs;
    ref_seqs = src.reference_sequences.dup;
//    sort!("a.length > b.length")(ref_seqs);
//    writefln("Refseqs: %s", ref_seqs);
    foreach (int refno, ReferenceSequenceInfo refseq; ref_seqs) {
        if ( canFind( refseq.name, "_" ) ) {
             continue;
         } 
//        stderr.writefln("Queuing task for contig %s", refseq.name);
        auto t = task!scanContig(_bamfilename, 
                                buildNormalizedPath(absolutePath(workdir)), 
                                refseq,
                                refno, 
                                cutoff_min, 
                                cutoff_max,
                                window_width, 
                                binsize,
                                h.q50, 
                                cast(uint) h.sd,
                                max_sd_degree,
                                n_threads);
        task_pool.put(t);
    }
    return 0;
}
