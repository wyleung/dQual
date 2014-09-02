module yamsvc.gzip;

import std.zlib;
import std.stdio;
import std.range;
import std.traits;

// source from: http://forum.dlang.org/thread/djhteyhpcnaskpabxijj@forum.dlang.org?page=2#post-cstqcyeajhyklvquhdwc:40forum.dlang.org

/*
import gzip;
import std.stdio;

void main() {
auto byLine = new GzipByLine("test.gz");
foreach(line; byLine)
   writeln(line);

auto gzipOutFile = new GzipOut("testout.gz");
gzipOutFile.compress("bla bla bla");
gzipOutFile.finish();
}
*/


class GzipInputRange {
  UnCompress uncompressObj;
  File f;
  auto CHUNKSIZE = 0x4000;
  ReturnType!(f.byChunk) chunkRange;
  bool exhausted;
  char[] uncompressedBuffer;
  size_t bufferIndex;

  this(string filename) {
    f = File(filename, "r");
    chunkRange = f.byChunk(CHUNKSIZE);
    uncompressObj = new UnCompress();
    load();
  }

  void load() {
    if(!chunkRange.empty) {
      auto raw = chunkRange.front.dup;
      chunkRange.popFront();
      uncompressedBuffer = cast(char[])uncompressObj.uncompress(raw);
      bufferIndex = 0;
    }
    else {
      if(!exhausted) {
        uncompressedBuffer = cast(char[])uncompressObj.flush();
        exhausted = true;
        bufferIndex = 0;
      }
      else
        uncompressedBuffer.length = 0;
    }
  }

  @property char front() {
    return uncompressedBuffer[bufferIndex];
  }

  void popFront() {
    bufferIndex += 1;
    if(bufferIndex >= uncompressedBuffer.length) {
      load();
      bufferIndex = 0;
    }
  }

  @property bool empty() {
    return uncompressedBuffer.length == 0;
  }
}

class GzipByLine {
  GzipInputRange range;
  char[] buf;

  this(string filename) {
    this.range = new GzipInputRange(filename);
    popFront();
  }

  @property bool empty() {
    return buf.length == 0;
  }

  void popFront() {
    buf.length = 0;
    while(!range.empty && range.front != '\n') {
      buf ~= range.front;
      range.popFront();
    }
    range.popFront();
  }

  string readln() {
  	return this.next();
  	}

  string next() {
  	this.popFront();
  	return this.front();
  	}

  string front() {
    return buf.idup;
  }
}

class GzipOut {
  Compress compressObj;
  File f;

  this(string filename) {
    f = File(filename, "w");
    compressObj = new Compress(HeaderFormat.gzip);
  }

  void compress(string s) {
    auto compressed = compressObj.compress(s.dup);
    f.rawWrite(compressed);
  }

  void finish() {
    auto compressed = compressObj.flush();
    f.rawWrite(compressed);
  }
}