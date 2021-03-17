/+dub.json:
{
	"name": "sir_stoch",
	"authors": [
		"Flávio Codeço Coelho"
	],
	"description": "Example usage of SIR model from epistochmodels",
	"copyright": "Copyright © 2020, Flávio Codeço Coelho",
	"license": "GPL-3",
    "dependencies": {
		"epistochmodels": "*"
	}
}
+/
import std.datetime.stopwatch;
import std.stdio;
import std.file;
import std.algorithm;
import std.range;
import core.exception : RangeError;
import std.typecons : tuple, Tuple;
import models;

/**
Save the simulation as a CSV file
*/
void save(string filename, string[] varnames, Tuple!(double[], int[][]) data)
{
    import std.array: join;
    import std.conv: text;
    File outf = File(filename, "w");
    auto head = join(varnames, ",");
    outf.writeln("t,", head);
    foreach (i, int[] row; data[1])
    {
        try
        {// writes t on column 0
            outf.writeln(text(data[0][i]), ",", row.map!text.join(","));
        }
        catch (RangeError e)
        {
            writefln("Some error %s: %s", e, row);
        }
    }
    outf.close();
}

void main(){
    double beta = 0.9;
    double gam = 0.1;
    int N = 100;
    int I0 = 10;
    double tf = 1000;
    auto sw = StopWatch(AutoStart.no);
    auto model = new SIR(N, beta, gam);
    model.initialize(N-I0, I0, 0);
    sw.start();
    auto sim = model.run(0, tf);
    sw.stop();
    save("sir_stoch.csv",["S","I","R"],sim);
    writefln("Time of the SIR run with N=%s: %s", N, sw.peek());
    writefln("Number of steps: %s", sim[0].length);
}
