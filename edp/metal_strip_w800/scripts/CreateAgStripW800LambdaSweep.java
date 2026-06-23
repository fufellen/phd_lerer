import com.comsol.model.Model;
import com.comsol.model.util.ModelUtil;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.io.PrintWriter;
import java.util.Locale;

public class CreateAgStripW800LambdaSweep {
    private static final int[] LAMBDAS_NM = new int[]{450, 500, 550, 600, 650, 700, 750, 800};
    private static final String[] EPS_AG = new String[]{
            "-6.078099+0.74594*i",
            "-8.491989+0.75842*i",
            "-11.126868+0.827824*i",
            "-13.904985+0.925288*i",
            "-17.021063+1.147584*i",
            "-20.42832+1.284248*i",
            "-23.969379+1.42042*i",
            "-27.953072+1.512654*i"
    };
    private static final String[] NEFF_GUESSES = new String[]{
            "2.39898973+0.205278903*i",
            "2.06858879+0.101873434*i",
            "1.88318726+0.0623954874*i",
            "1.77102486+0.0432107061*i",
            "1.69035403+0.034190244*i",
            "1.63321641+0.0252388481*i",
            "1.59336098+0.0192723921*i",
            "1.56145898+0.014141931*i"
    };

    private static final String FIELD_TOTAL =
            "intopAg(ewfd.normE^2)+intopSub(ewfd.normE^2)+intopAirL(ewfd.normE^2)+intopAirR(ewfd.normE^2)+intopAirTop(ewfd.normE^2)";

    public static void main(String[] args) throws IOException {
        writeSummary("comsol_w800_lambda_sweep_status.txt",
                "Starting COMSOL wavelength sweep for Ag strip W=800 nm.", false);
        for (int i = 0; i < LAMBDAS_NM.length; i++) {
            runLambda(LAMBDAS_NM[i], EPS_AG[i], NEFF_GUESSES[i]);
        }
        writeSummary("comsol_w800_lambda_sweep_status.txt",
                "Completed COMSOL wavelength sweep for Ag strip W=800 nm.", true);
    }

    public static Model runLambda(int lambdaNm, String epsAg, String neffGuess) throws IOException {
        String label = String.format(Locale.US, "lam%04dnm", lambdaNm);
        String modelTag = "Model_" + label;
        double lambdaUm = lambdaNm / 1000.0;
        String lambdaUmText = String.format(Locale.US, "%.6f[um]", lambdaUm);

        writeSummary("comsol_w800_lambda_sweep_status.txt",
                "Running lambda " + lambdaNm + " nm with eps_Ag=" + epsAg + " and shift " + neffGuess + ".", true);

        Model model = ModelUtil.create(modelTag);
        model.modelPath(System.getProperty("user.dir"));
        model.label("ag_strip_w0800_" + label + "_mode_2d.mph");

        model.param().set("lambda0", lambdaUmText);
        model.param().set("f0", "c_const/lambda0");
        model.param().set("n_sub", "1.45");
        model.param().set("eps_sub", "n_sub^2");
        model.param().set("eps_air", "1");
        model.param().set("eps_Ag", epsAg);
        model.param().set("t_Ag", "0.020[um]");
        model.param().set("w_Ag", "0.800[um]");
        model.param().set("w_box", "6.0[um]");
        model.param().set("h_air", "2.0[um]");
        model.param().set("h_sub", "2.0[um]");
        model.param().set("neff_guess", neffGuess);

        model.component().create("comp1", true);
        model.component("comp1").geom().create("geom1", 2);
        model.component("comp1").geom("geom1").lengthUnit("um");

        addRectangle(model, "rSub", "Glass substrate, n=1.45",
                new String[]{"-w_box/2", "-h_sub"}, new String[]{"w_box", "h_sub"});
        addRectangle(model, "rAirL", "Air, left of strip",
                new String[]{"-w_box/2", "0"}, new String[]{"(w_box-w_Ag)/2", "h_air"});
        addRectangle(model, "rAirR", "Air, right of strip",
                new String[]{"w_Ag/2", "0"}, new String[]{"(w_box-w_Ag)/2", "h_air"});
        addRectangle(model, "rAirTop", "Air over strip",
                new String[]{"-w_Ag/2", "t_Ag"}, new String[]{"w_Ag", "h_air-t_Ag"});
        addRectangle(model, "rAg", "Ag strip, 20 nm x 800 nm",
                new String[]{"-w_Ag/2", "0"}, new String[]{"w_Ag", "t_Ag"});
        model.component("comp1").geom("geom1").run();

        addMaterial(model, "matSub", "Glass substrate, n=1.45", "geom1_rSub_dom", "eps_sub");
        addMaterial(model, "matAg", "Ag from Ag_.c at " + lambdaNm + " nm", "geom1_rAg_dom", "eps_Ag");
        addMaterial(model, "matAirL", "Air, left", "geom1_rAirL_dom", "eps_air");
        addMaterial(model, "matAirR", "Air, right", "geom1_rAirR_dom", "eps_air");
        addMaterial(model, "matAirTop", "Air, top", "geom1_rAirTop_dom", "eps_air");

        addIntegration(model, "intopAg", "geom1_rAg_dom");
        addIntegration(model, "intopSub", "geom1_rSub_dom");
        addIntegration(model, "intopAirL", "geom1_rAirL_dom");
        addIntegration(model, "intopAirR", "geom1_rAirR_dom");
        addIntegration(model, "intopAirTop", "geom1_rAirTop_dom");

        model.component("comp1").physics().create("ewfd", "ElectromagneticWavesFrequencyDomain", "geom1");
        createMesh(model);
        createStudy(model, label);
        model.study("std1").run();

        exportModeTable(model, label);
        model.save("ag_strip_w0800_" + label + "_mode_2d.mph");

        writeSummary("comsol_w800_lambda_sweep_status.txt",
                "Completed lambda " + lambdaNm + " nm.", true);
        ModelUtil.remove(modelTag);
        return model;
    }

    private static void addRectangle(Model model, String tag, String label, String[] pos, String[] size) {
        model.component("comp1").geom("geom1").create(tag, "Rectangle");
        model.component("comp1").geom("geom1").feature(tag).label(label);
        model.component("comp1").geom("geom1").feature(tag).set("pos", pos);
        model.component("comp1").geom("geom1").feature(tag).set("size", size);
        model.component("comp1").geom("geom1").feature(tag).set("selresult", "on");
        model.component("comp1").geom("geom1").feature(tag).set("selresultshow", "dom");
    }

    private static void addMaterial(Model model, String tag, String label, String selection, String eps) {
        model.component("comp1").material().create(tag, "Common");
        model.component("comp1").material(tag).label(label);
        model.component("comp1").material(tag).selection().named(selection);
        model.component("comp1").material(tag).propertyGroup("def").set("relpermittivity",
                new String[]{eps, "0", "0", "0", eps, "0", "0", "0", eps});
    }

    private static void addIntegration(Model model, String tag, String selection) {
        model.component("comp1").cpl().create(tag, "Integration");
        model.component("comp1").cpl(tag).selection().named(selection);
    }

    private static void createMesh(Model model) {
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
    }

    private static void createStudy(Model model, String label) {
        model.study().create("std1");
        model.study("std1").label("2D Mode Analysis near EDP estimate, " + label);
        model.study("std1").create("mode", "ModeAnalysis");
        model.study("std1").feature("mode").set("modeFreq", "f0");
        model.study("std1").feature("mode").set("neigsactive", true);
        model.study("std1").feature("mode").set("neigs", 20);
        model.study("std1").feature("mode").set("shiftactive", true);
        model.study("std1").feature("mode").set("shift", "neff_guess");
    }

    private static void exportModeTable(Model model, String label) {
        String airExpr = "intopAirL(ewfd.normE^2)+intopAirR(ewfd.normE^2)+intopAirTop(ewfd.normE^2)";
        model.result().table().create("tbl1", "Table");
        model.result().table("tbl1").label("Ag strip W=800 nm modes at " + label);
        model.result().numerical().create("gev1", "EvalGlobal");
        model.result().numerical("gev1").label("Mode quantities and field localization proxies");
        model.result().numerical("gev1").set("expr",
                new String[]{
                        "ewfd.neff",
                        "real(ewfd.neff)",
                        "imag(ewfd.neff)",
                        "imag(ewfd.neff)/real(ewfd.neff)",
                        "ewfd.beta",
                        "ewfd.dampzdB",
                        "intopAg(ewfd.normE^2)",
                        "intopSub(ewfd.normE^2)",
                        airExpr,
                        "intopAg(ewfd.normE^2)/(" + FIELD_TOTAL + ")",
                        "intopSub(ewfd.normE^2)/(" + FIELD_TOTAL + ")",
                        "(" + airExpr + ")/(" + FIELD_TOTAL + ")"
                });
        model.result().numerical("gev1").set("unit",
                new String[]{"1", "1", "1", "1", "rad/um", "dB/um", "V^2", "V^2", "V^2", "1", "1", "1"});
        model.result().numerical("gev1").set("descr",
                new String[]{
                        "Effective mode index",
                        "Re(neff)",
                        "Im(neff)",
                        "Im/Re",
                        "Propagation constant beta",
                        "Attenuation constant",
                        "Integral |E|^2 in Ag",
                        "Integral |E|^2 in substrate",
                        "Integral |E|^2 in air",
                        "Ag field fraction proxy",
                        "Substrate field fraction proxy",
                        "Air field fraction proxy"
                });
        model.result().numerical("gev1").set("table", "tbl1");
        model.result().numerical("gev1").setResult();
        model.result().table("tbl1").save("comsol_ag_strip_w800_lambda_sweep_" + label + "_modes.csv");
    }

    private static void writeSummary(String file, String text, boolean append) throws IOException {
        PrintWriter out = new PrintWriter(new OutputStreamWriter(new FileOutputStream(file, append), "UTF-8"));
        try {
            out.println(text);
        } finally {
            out.close();
        }
    }
}
