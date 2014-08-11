module yamsvc.utils;

import std.algorithm : count, max, reduce;
import std.math : abs;
import bio.bam.read;

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

struct BedRecord {
	string chromosome;
	ulong start;
	ulong end;
	string infofields;
}

class DiscordantRegion {
	ulong start = 0;
	ulong end = 0;
	short chromosome=-1;
	ushort[] coverages_discordant;
	ushort[] coverages;
	bool working;
	
	void addCoverage( ushort cov ) {
		this.coverages ~= cov;
	}
	
	void addCoverageDiscordant( ushort cov ) {
		this.coverages_discordant ~= cov;
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
	
	void resetCoverage() {
		/*
			Empty the coverages for a new region run
		*/
		this.coverages = [];
		this.coverages_discordant = [];
	}
	
	ulong length() {
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
}
