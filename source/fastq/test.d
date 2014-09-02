#!/usr/bin/env rdmd

/*
 @author: Wai Yi Leung
 @copyright: 2013-2014 Wai Yi Leung
 @contact: w.y.leung@e-sensei.nl
 @license: MIT Licence
*/

import yamsvc.gzip;
import std.stdio;
import std.file;
import std.process;

int main(string[] args) {
	
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

    string line;    
    while ( (line = file.readln()) !is null ) {
        string name = line;
        string sequence = file.readln();
        string holder = file.readln();
        string quality = file.readln();
        writefln("Name: %s", name);
    }
	
	
//	auto byLine = new GzipByLine(args[1]);
//	foreach(line; byLine) {
//		writeln(line);
//		writeln(byLine.next());
//		writeln(byLine.next());
//		writeln(byLine.next());
//		writeln("Next!");
//	}
	return 0;
}