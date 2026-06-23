import com.comsol.model.Model;
import com.comsol.model.util.ModelUtil;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.io.PrintWriter;

public class CreateAgStripW800Mode2D {
    public static Model run() throws IOException {
        Model model = ModelUtil.create("Model");
        model.modelPath(System.getProperty("user.dir"));
        model.label("ag_strip_w800_500nm_mode_2d.mph");

        model.param().set("lambda0", "0.50[um]");
        model.param().set("f0", "c_const/lambda0");
        model.param().set("n_sub", "1.45");
        model.param().set("eps_sub", "n_sub^2");
        model.param().set("eps_air", "1");
        model.param().set("eps_Ag", "-9.4079+0.7982*i");
        model.param().set("t_Ag", "0.020[um]");
        model.param().set("w_Ag", "0.800[um]");
        model.param().set("w_box", "6.0[um]");
        model.param().set("h_air", "2.0[um]");
        model.param().set("h_sub", "2.0[um]");
        model.param().set("neff_guess", "2.0686+0.1018*i");

        model.component().create("comp1", true);
        model.component("comp1").geom().create("geom1", 2);
        model.component("comp1").geom("geom1").lengthUnit("um");

        model.component("comp1").geom("geom1").create("rSub", "Rectangle");
        model.component("comp1").geom("geom1").feature("rSub").label("Glass substrate, n=1.45");
        model.component("comp1").geom("geom1").feature("rSub").set("pos", new String[]{"-w_box/2", "-h_sub"});
        model.component("comp1").geom("geom1").feature("rSub").set("size", new String[]{"w_box", "h_sub"});
        model.component("comp1").geom("geom1").feature("rSub").set("selresult", "on");
        model.component("comp1").geom("geom1").feature("rSub").set("selresultshow", "dom");

        model.component("comp1").geom("geom1").create("rAirL", "Rectangle");
        model.component("comp1").geom("geom1").feature("rAirL").label("Air, left of strip");
        model.component("comp1").geom("geom1").feature("rAirL").set("pos", new String[]{"-w_box/2", "0"});
        model.component("comp1").geom("geom1").feature("rAirL").set("size", new String[]{"(w_box-w_Ag)/2", "h_air"});
        model.component("comp1").geom("geom1").feature("rAirL").set("selresult", "on");
        model.component("comp1").geom("geom1").feature("rAirL").set("selresultshow", "dom");

        model.component("comp1").geom("geom1").create("rAirR", "Rectangle");
        model.component("comp1").geom("geom1").feature("rAirR").label("Air, right of strip");
        model.component("comp1").geom("geom1").feature("rAirR").set("pos", new String[]{"w_Ag/2", "0"});
        model.component("comp1").geom("geom1").feature("rAirR").set("size", new String[]{"(w_box-w_Ag)/2", "h_air"});
        model.component("comp1").geom("geom1").feature("rAirR").set("selresult", "on");
        model.component("comp1").geom("geom1").feature("rAirR").set("selresultshow", "dom");

        model.component("comp1").geom("geom1").create("rAirTop", "Rectangle");
        model.component("comp1").geom("geom1").feature("rAirTop").label("Air over strip");
        model.component("comp1").geom("geom1").feature("rAirTop").set("pos", new String[]{"-w_Ag/2", "t_Ag"});
        model.component("comp1").geom("geom1").feature("rAirTop").set("size", new String[]{"w_Ag", "h_air-t_Ag"});
        model.component("comp1").geom("geom1").feature("rAirTop").set("selresult", "on");
        model.component("comp1").geom("geom1").feature("rAirTop").set("selresultshow", "dom");

        model.component("comp1").geom("geom1").create("rAg", "Rectangle");
        model.component("comp1").geom("geom1").feature("rAg").label("Ag strip, 20 nm x 800 nm");
        model.component("comp1").geom("geom1").feature("rAg").set("pos", new String[]{"-w_Ag/2", "0"});
        model.component("comp1").geom("geom1").feature("rAg").set("size", new String[]{"w_Ag", "t_Ag"});
        model.component("comp1").geom("geom1").feature("rAg").set("selresult", "on");
        model.component("comp1").geom("geom1").feature("rAg").set("selresultshow", "dom");
        model.component("comp1").geom("geom1").run();

        model.component("comp1").material().create("matSub", "Common");
        model.component("comp1").material("matSub").label("Glass substrate, n=1.45");
        model.component("comp1").material("matSub").selection().named("geom1_rSub_dom");
        model.component("comp1").material("matSub").propertyGroup("def").set("relpermittivity",
                new String[]{"eps_sub", "0", "0", "0", "eps_sub", "0", "0", "0", "eps_sub"});

        model.component("comp1").material().create("matAg", "Common");
        model.component("comp1").material("matAg").label("Ag from Ag_.c at 500 nm");
        model.component("comp1").material("matAg").selection().named("geom1_rAg_dom");
        model.component("comp1").material("matAg").propertyGroup("def").set("relpermittivity",
                new String[]{"eps_Ag", "0", "0", "0", "eps_Ag", "0", "0", "0", "eps_Ag"});

        model.component("comp1").material().create("matAirL", "Common");
        model.component("comp1").material("matAirL").label("Air, left");
        model.component("comp1").material("matAirL").selection().named("geom1_rAirL_dom");
        model.component("comp1").material("matAirL").propertyGroup("def").set("relpermittivity",
                new String[]{"eps_air", "0", "0", "0", "eps_air", "0", "0", "0", "eps_air"});

        model.component("comp1").material().create("matAirR", "Common");
        model.component("comp1").material("matAirR").label("Air, right");
        model.component("comp1").material("matAirR").selection().named("geom1_rAirR_dom");
        model.component("comp1").material("matAirR").propertyGroup("def").set("relpermittivity",
                new String[]{"eps_air", "0", "0", "0", "eps_air", "0", "0", "0", "eps_air"});

        model.component("comp1").material().create("matAirTop", "Common");
        model.component("comp1").material("matAirTop").label("Air, top");
        model.component("comp1").material("matAirTop").selection().named("geom1_rAirTop_dom");
        model.component("comp1").material("matAirTop").propertyGroup("def").set("relpermittivity",
                new String[]{"eps_air", "0", "0", "0", "eps_air", "0", "0", "0", "eps_air"});

        model.component("comp1").physics().create("ewfd", "ElectromagneticWavesFrequencyDomain", "geom1");

        model.component("comp1").mesh().create("mesh1");
        model.component("comp1").mesh("mesh1").feature("size").set("custom", "on");
        model.component("comp1").mesh("mesh1").feature("size").set("hmax", "0.06");
        model.component("comp1").mesh("mesh1").feature("size").set("hmin", "0.001");
        model.component("comp1").mesh("mesh1").feature("size").set("hgrad", "1.2");

        model.component("comp1").mesh("mesh1").feature().create("sizeAg", "Size");
        model.component("comp1").mesh("mesh1").feature("sizeAg").label("Fine mesh in Ag strip");
        model.component("comp1").mesh("mesh1").feature("sizeAg").selection().named("geom1_rAg_dom");
        model.component("comp1").mesh("mesh1").feature("sizeAg").set("custom", "on");
        model.component("comp1").mesh("mesh1").feature("sizeAg").set("hmax", "0.003");
        model.component("comp1").mesh("mesh1").feature("sizeAg").set("hmin", "0.0005");
        model.component("comp1").mesh("mesh1").feature("sizeAg").set("hgrad", "1.15");

        model.component("comp1").mesh("mesh1").feature().create("ftri1", "FreeTri");
        model.component("comp1").mesh("mesh1").run();

        model.study().create("std1");
        model.study("std1").label("2D Mode Analysis near EDP estimate");
        model.study("std1").create("mode", "ModeAnalysis");
        model.study("std1").feature("mode").set("modeFreq", "f0");
        model.study("std1").feature("mode").set("neigsactive", true);
        model.study("std1").feature("mode").set("neigs", 20);
        model.study("std1").feature("mode").set("shiftactive", true);
        model.study("std1").feature("mode").set("shift", "neff_guess");

        model.study("std1").run();

        model.result().table().create("tbl1", "Table");
        model.result().table("tbl1").label("Ag strip W=800 nm modes at 500 nm");
        model.result().numerical().create("gev1", "EvalGlobal");
        model.result().numerical("gev1").label("Mode quantities");
        model.result().numerical("gev1").set("expr",
                new String[]{"ewfd.neff", "real(ewfd.neff)", "imag(ewfd.neff)", "imag(ewfd.neff)/real(ewfd.neff)", "ewfd.beta", "ewfd.dampzdB"});
        model.result().numerical("gev1").set("unit",
                new String[]{"1", "1", "1", "1", "rad/um", "dB/um"});
        model.result().numerical("gev1").set("descr",
                new String[]{"Effective mode index", "Re(neff)", "Im(neff)", "Im/Re", "Propagation constant beta", "Attenuation constant"});
        model.result().numerical("gev1").set("table", "tbl1");
        model.result().numerical("gev1").setResult();
        model.result().table("tbl1").save("comsol_ag_strip_w800_500nm_modes.csv");

        model.result().create("pg1", "PlotGroup2D");
        model.result("pg1").label("Ag strip W=800 nm: electric field norm");
        model.result("pg1").feature().create("surf1", "Surface");
        model.result("pg1").feature("surf1").set("expr", "ewfd.normE");
        model.result("pg1").feature("surf1").set("unit", "V/m");

        model.save("ag_strip_w800_500nm_mode_2d.mph");
        writeSummary("comsol_ag_strip_w800_500nm_status.txt",
                "COMSOL 2D Mode Analysis completed for Ag strip W=800 nm, t=20 nm, substrate n=1.45, lambda=500 nm.");
        return model;
    }

    private static void writeSummary(String file, String text) throws IOException {
        PrintWriter out = new PrintWriter(new OutputStreamWriter(new FileOutputStream(file), "UTF-8"));
        try {
            out.println(text);
        } finally {
            out.close();
        }
    }

    public static void main(String[] args) throws IOException {
        run();
    }
}
