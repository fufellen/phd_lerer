import com.comsol.model.Model;
import com.comsol.model.util.ModelUtil;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.io.PrintWriter;
import java.util.Locale;

/**
 * The butt joint between the silicon feed and the Sb2S3/Au plasmonic phase shifter, driven from
 * EITHER end.  The geometry is the one of RunMatching.java (J_junction): half structure closed by
 * a magnetic wall at x = 0, the 10 nm gold film as a transition boundary condition with the
 * calibrated 184.356 nm ridge, numeric ports on both ends, scattering boundaries elsewhere.
 *
 * What is new here:
 *   - MR_DIR = fwd excites port 1 (the feed), MR_DIR = rev excites port 2 (the phase shifter),
 *     with the excitation of the OTHER port switched off explicitly;
 *   - the S-parameters of the numeric ports are written for every run, so the modal
 *     transmission of the joint follows from |S21| or |S12| and the decay of the port mode
 *     over the uniform section, without any fit to a power profile;
 *   - the forward power is sampled on planes on BOTH sides of the joint, so the reverse run
 *     shows how the power arriving in the lossless feed settles;
 *   - the port mode indices are written too, so the decay used in the extraction is the model's
 *     own and not borrowed from a 2D section.
 *
 * Schemes (MR_ONLY): abrupt (736 x 240 nm feed on SiO2, 120 nm step), shared450, sync405.
 * MR_PCM: 2.712 amorphous, 3.308 crystalline.  MR_MESH scales every mesh size.  MR_LFEED and
 * MR_LPCM override the section lengths (defaults 0.9 and 3.0 um, as in RunMatching).
 * MR_IMAGE=1 exports the side view and two top views, MR_SAVE=1 keeps the solved model.
 */
public class RunReverse {

    private static final String STATUS = "reverse_status.txt";
    private static PrintWriter out;
    private static PrintWriter prof;
    private static double meshScale = 1.0;
    private static double lFeed = 0.9;
    private static double lPcm = 3.0;

    public static void main(String[] args) throws IOException {
        write(STATUS, "START reverse", false);
        boolean fresh = !new java.io.File("reverse_results.csv").exists();
        out = utf8("reverse_results.csv", !fresh);
        if (fresh) {
            out.println("case,scheme,dir,n_pcm,mesh,l_feed_um,l_pcm_um,elements,"
                    + "neff1_re,neff1_im,neff2_re,neff2_im,"
                    + "S11_re,S11_im,S21_re,S21_im,S12_re,S12_im,S22_re,S22_im,"
                    + "p_feed_end,p_pcm_end");
        }
        out.flush();
        boolean freshP = !new java.io.File("reverse_profiles.csv").exists();
        prof = utf8("reverse_profiles.csv", !freshP);
        if (freshP) {
            prof.println("case,z_um,p_z");
        }
        prof.flush();

        if (System.getenv("MR_MESH") != null) {
            meshScale = Double.parseDouble(System.getenv("MR_MESH"));
        }
        lFeed = env("MR_LFEED", 0.9);
        lPcm = env("MR_LPCM", 3.0);
        String only = System.getenv("MR_ONLY");
        if (only == null) {
            only = "sync405";
        }
        String dir = System.getenv("MR_DIR");
        if (dir == null) {
            dir = "fwd";
        }
        boolean forward = dir.equals("fwd");
        double nPcm = env("MR_PCM", 2.712);
        boolean cryst = nPcm > 3.0;

        // shifts for the port mode searches: the feed index depends on the scheme and, on the
        // shared buffer, on the PCM state; the phase-shifter index on the state only
        String shift2 = cryst ? "2.2772+1.50e-2*i" : "1.8877+2.46e-2*i";
        if (only.equals("abrupt")) {
            run("abrupt", forward, nPcm, false, 0.736, 1.8877, shift2);
        } else if (only.equals("shared450")) {
            run("shared450", forward, nPcm, true, 0.450, cryst ? env("MR_NFEED", 2.23) : 1.9326,
                    shift2);
        } else if (only.equals("sync405")) {
            run("sync405", forward, nPcm, true, 0.405, cryst ? env("MR_NFEED", 2.175) : 1.8877,
                    shift2);
        } else {
            write(STATUS, "unknown scheme " + only, true);
        }
        out.close();
        prof.close();
        write(STATUS, "DONE", true);
    }

    private static double env(String name, double dflt) {
        String v = System.getenv(name);
        return v == null ? dflt : Double.parseDouble(v);
    }

    private static void run(String scheme, boolean forward, double nPcm, boolean shared,
                            double wFeed, double shift1, String shift2) throws IOException {
        String caseName = scheme + "_" + (forward ? "fwd" : "rev") + (nPcm > 3.0 ? "_c" : "_a")
                + (meshScale != 1.0 ? "_m" + Math.round(meshScale * 100) : "")
                + (lFeed != 0.9 ? "_lf" + Math.round(lFeed * 100) : "");
        String tag = "MR_" + scheme;
        write(STATUS, "CASE " + caseName, true);
        Model model = ModelUtil.create(tag);
        model.modelPath(System.getProperty("user.dir"));

        double wHalf = 1.0;
        double tPcm = 0.120, tLoad = 0.184356, tWg = 0.240, wCore = 0.450;

        model.param().set("lambda0", "1.55[um]");
        model.param().set("f0", "c_const/lambda0");
        model.param().set("n_sio2", "1.444");
        model.param().set("n_si", "3.478");
        // loss is a NEGATIVE imaginary permittivity in the frequency domain
        model.param().set("eps_au", "(0.6389-11.1748*i)^2");
        model.param().set("n_pcm", String.format(Locale.US, "%.4f", nPcm));
        model.param().set("t_au", "0.010[um]");

        model.component().create("comp1", true);
        model.component("comp1").geom().create("geom1", 3);
        model.component("comp1").geom("geom1").lengthUnit("um");
        double zPcm0 = shared ? -lFeed : 0.0;
        blk(model, "bSub", 0, -tPcm - 0.70, -lFeed, wHalf, 0.70, lFeed + lPcm);
        blk(model, "bAir", 0, -tPcm, -lFeed, wHalf, 0.70 + tPcm, lFeed + lPcm);
        blk(model, "bPcm", 0, -tPcm, zPcm0, 0.400, tPcm, lPcm - zPcm0);
        if (shared) {
            blk(model, "bWg", 0, 0.0, -lFeed, wFeed / 2.0, tLoad, lFeed);
        } else {
            blk(model, "bWg", 0, -tPcm, -lFeed, wFeed / 2.0, tWg, lFeed);
        }
        blk(model, "bCore0", 0, 0.0, 0.0, wCore / 2.0, tLoad, lPcm);
        model.component("comp1").geom("geom1").run();

        matN(model, "mAir", "geom1_bAir_dom", "1");
        matN(model, "mSub", "geom1_bSub_dom", "n_sio2");
        matN(model, "mPcm", "geom1_bPcm_dom", "n_pcm");
        matN(model, "mWg", "geom1_bWg_dom", "n_si");
        matN(model, "mCore0", "geom1_bCore0_dom", "n_si");

        model.component("comp1").physics().create("ewfd",
                "ElectromagneticWavesFrequencyDomain", "geom1");

        boxSel(model, "selSym", 2, -1e-3, 1e-3, -1e4, 1e4, -1e4, 1e4);
        model.component("comp1").physics("ewfd").create("pmc1", "PerfectMagneticConductor", 2);
        model.component("comp1").physics("ewfd").feature("pmc1").selection().named("selSym");

        boxSel(model, "selOutX", 2, wHalf - 1e-3, wHalf + 1e-3, -1e4, 1e4, -1e4, 1e4);
        boxSel(model, "selOutYlo", 2, -1e4, 1e4, -tPcm - 0.70 - 1e-3, -tPcm - 0.70 + 1e-3,
                -1e4, 1e4);
        boxSel(model, "selOutYhi", 2, -1e4, 1e4, 0.70 - 1e-3, 0.70 + 1e-3, -1e4, 1e4);
        model.component("comp1").physics("ewfd").create("sbc1", "Scattering", 2);
        model.component("comp1").physics("ewfd").feature("sbc1").selection().set(
                union(model, "selOutX", "selOutYlo", "selOutYhi"));

        // the gold: the y = 0 face under the core, from the joint to the far port
        boxSel(model, "selAu0", 2, -1e-3, wCore / 2.0 + 1e-3, -1e-3, 1e-3, -1e-4, lPcm + 1e-4);
        int[] auFaces = model.component("comp1").selection("selAu0").entities(2);
        model.component("comp1").selection().create("selAu", "Explicit");
        model.component("comp1").selection("selAu").geom(2);
        model.component("comp1").selection("selAu").set(auFaces);
        write(STATUS, "  gold faces: " + auFaces.length, true);
        model.component("comp1").physics("ewfd").create("tbc1", "TransitionBoundaryCondition", 2);
        model.component("comp1").physics("ewfd").feature("tbc1").selection().named("selAu");
        model.component("comp1").physics("ewfd").feature("tbc1").set("d", "t_au");
        model.component("comp1").material().create("mAuSheet", "Common");
        model.component("comp1").material("mAuSheet").label("Gold sheet");
        model.component("comp1").material("mAuSheet").selection().geom("geom1", 2);
        model.component("comp1").material("mAuSheet").selection().named("selAu");
        model.component("comp1").material("mAuSheet").propertyGroup("def")
                .set("relpermittivity", new String[]{"eps_au"});
        model.component("comp1").material("mAuSheet").propertyGroup("def")
                .set("relpermeability", new String[]{"1"});
        model.component("comp1").material("mAuSheet").propertyGroup("def")
                .set("electricconductivity", new String[]{"0"});

        boxSel(model, "selP1", 2, -1e4, 1e4, -1e4, 1e4, -lFeed - 1e-3, -lFeed + 1e-3);
        boxSel(model, "selP2", 2, -1e4, 1e4, -1e4, 1e4, lPcm - 1e-3, lPcm + 1e-3);
        port(model, "port1", "selP1", forward, "1");
        port(model, "port2", "selP2", !forward, "2");

        model.component("comp1").mesh().create("mesh1");
        model.component("comp1").mesh("mesh1").feature("size").set("custom", "on");
        model.component("comp1").mesh("mesh1").feature("size").set("hmax",
                String.valueOf(0.16 * meshScale));
        model.component("comp1").mesh("mesh1").feature("size").set("hmin",
                String.valueOf(0.006 * meshScale));
        model.component("comp1").mesh("mesh1").feature("size").set("hgrad", "1.5");
        sizeDom(model, "szWg", "geom1_bWg_dom", 0.070 * meshScale);
        sizeDom(model, "szPcm", "geom1_bPcm_dom", 0.060 * meshScale);
        sizeDom(model, "szCore0", "geom1_bCore0_dom", 0.060 * meshScale);
        sizeBnd(model, "szAu", "selAu", 0.035 * meshScale);
        model.component("comp1").mesh("mesh1").feature().create("ftet1", "FreeTet");
        model.component("comp1").mesh("mesh1").run();
        int nel = model.component("comp1").mesh("mesh1").getNumElem();
        write(STATUS, "  meshed, elements=" + nel, true);

        model.study().create("std1");
        model.study("std1").create("bma1", "BoundaryModeAnalysis");
        model.study("std1").feature("bma1").set("PortName", "1");
        model.study("std1").feature("bma1").set("modeFreq", "f0");
        model.study("std1").feature("bma1").set("neigsactive", true);
        model.study("std1").feature("bma1").set("neigs", 1);
        model.study("std1").feature("bma1").set("shiftactive", true);
        model.study("std1").feature("bma1").set("shift", String.format(Locale.US, "%.4f", shift1));
        model.study("std1").create("bma2", "BoundaryModeAnalysis");
        model.study("std1").feature("bma2").set("PortName", "2");
        model.study("std1").feature("bma2").set("modeFreq", "f0");
        model.study("std1").feature("bma2").set("neigsactive", true);
        model.study("std1").feature("bma2").set("neigs", 1);
        model.study("std1").feature("bma2").set("shift", shift2);
        model.study("std1").feature("bma2").set("shiftactive", true);
        model.study("std1").create("freq", "Frequency");
        model.study("std1").feature("freq").set("plist", "f0");
        model.study("std1").run();
        write(STATUS, "  solved", true);

        // port mode indices and S-parameters, each in its own evaluation so that a name that
        // does not evaluate costs one NaN and not the whole row
        double[] neff1 = global2(model, "ewfd.neff_1");
        double[] neff2 = global2(model, "ewfd.neff_2");
        if (Double.isNaN(neff1[0])) {
            double[] b1 = global2(model, "ewfd.beta_1");
            double[] b2 = global2(model, "ewfd.beta_2");
            double k0 = 2.0 * Math.PI / 1.55e-6;
            neff1 = new double[]{b1[0] / k0, b1[1] / k0};
            neff2 = new double[]{b2[0] / k0, b2[1] / k0};
        }
        double[] s11 = global2(model, "ewfd.S11");
        double[] s21 = global2(model, "ewfd.S21");
        double[] s12 = global2(model, "ewfd.S12");
        double[] s22 = global2(model, "ewfd.S22");
        write(STATUS, String.format(Locale.US,
                "  neff1=%.6f%+.3ei neff2=%.6f%+.3ei S11=%.4e%+.4ei S21=%.4e%+.4ei S12=%.4e%+.4ei S22=%.4e%+.4ei",
                neff1[0], neff1[1], neff2[0], neff2[1], s11[0], s11[1], s21[0], s21[1],
                s12[0], s12[1], s22[0], s22[1]), true);

        // forward power on planes on both sides of the joint; negative in the reverse run
        java.util.List<Double> zl = new java.util.ArrayList<Double>();
        for (double z = -lFeed + 0.1; z < -0.15; z += (lFeed > 1.2 ? 0.3 : 0.25)) {
            zl.add(Math.round(z * 1000.0) / 1000.0);
        }
        zl.add(-0.10);
        zl.add(0.15);
        for (double z = 0.5; z < lPcm - 0.15; z += 0.5) {
            zl.add(z);
        }
        zl.add(lPcm - 0.2);
        StringBuilder sb = new StringBuilder("  ");
        double pFirst = Double.NaN, pLast = Double.NaN;
        int k = 0;
        for (double z : zl) {
            double p = power(model, "cp" + (k++), z);
            if (Double.isNaN(pFirst)) {
                pFirst = p;
            }
            pLast = p;
            prof.println(String.format(Locale.US, "%s,%.4f,%.9e", caseName, z, p));
            sb.append(String.format(Locale.US, "P(%.2f)=%.4f ", z, p));
        }
        prof.flush();
        write(STATUS, sb.toString(), true);
        out.println(String.format(Locale.US,
                "%s,%s,%s,%.4f,%.2f,%.3f,%.3f,%d,%.7f,%.4e,%.7f,%.4e,"
                + "%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.9e,%.9e",
                caseName, scheme, forward ? "fwd" : "rev", nPcm, meshScale, lFeed, lPcm, nel,
                neff1[0], neff1[1], neff2[0], neff2[1],
                s11[0], s11[1], s21[0], s21[1], s12[0], s12[1], s22[0], s22[1], pFirst, pLast));
        out.flush();
        if (System.getenv("MR_IMAGE") != null) {
            try {
                view(model, caseName, "side", "zy", "quickx", 0.0);
                view(model, caseName, "top", "zx", "quicky", 0.005);
                view(model, caseName, "buf", "zx", "quicky", -0.060);
            } catch (Exception ex) {
                write(STATUS, "  image failed: " + ex, true);
            }
        }
        if (System.getenv("MR_SAVE") != null) {
            try {
                model.save(System.getProperty("user.dir") + java.io.File.separator
                        + "reverse_" + caseName + ".mph");
                write(STATUS, "  saved reverse_" + caseName + ".mph", true);
            } catch (Exception ex) {
                write(STATUS, "  save failed: " + ex, true);
            }
        }
        write(STATUS, "CASE " + caseName + " OK", true);
        ModelUtil.remove(tag);
    }

    /** Real and imaginary part of a global expression, NaN when it does not evaluate. */
    private static double[] global2(Model model, String expr) {
        String num = "gl_" + expr.replace('.', '_');
        try {
            model.result().numerical().create(num, "EvalGlobal");
            model.result().numerical(num).set("expr", new String[]{expr});
            double re = model.result().numerical(num).getReal()[0][0];
            double im = model.result().numerical(num).getImag()[0][0];
            model.result().numerical().remove(num);
            return new double[]{re, im};
        } catch (Exception ex) {
            try {
                model.result().numerical().remove(num);
            } catch (Exception ignored) {
                // nothing to remove
            }
            try {
                write(STATUS, "  " + expr + " does not evaluate: " + ex.getMessage(), true);
            } catch (IOException ignored) {
                // status file unavailable
            }
            return new double[]{Double.NaN, Double.NaN};
        }
    }

    private static double power(Model model, String tag, double z) {
        model.result().dataset().create(tag, "CutPlane");
        model.result().dataset(tag).set("quickplane", "xy");
        model.result().dataset(tag).set("quickz", String.valueOf(z));
        String num = "si_" + tag;
        model.result().numerical().create(num, "IntSurface");
        model.result().numerical(num).set("data", tag);
        model.result().numerical(num).set("expr", new String[]{"ewfd.Poavz"});
        double v = model.result().numerical(num).getReal()[0][0];
        model.result().numerical().remove(num);
        model.result().dataset().remove(tag);
        return v;
    }

    private static void view(Model model, String scheme, String name, String plane,
                             String offsetProp, double offset) {
        model.result().dataset().create("cpImg", "CutPlane");
        model.result().dataset("cpImg").set("quickplane", plane);
        model.result().dataset("cpImg").set(offsetProp, String.valueOf(offset));
        model.result().create("pgImg", "PlotGroup2D");
        model.result("pgImg").set("data", "cpImg");
        model.result("pgImg").create("surf1", "Surface");
        model.result("pgImg").feature("surf1").set("expr", "ewfd.normE");
        model.result("pgImg").feature("surf1").set("colortable", "Inferno");
        model.result("pgImg").run();
        model.result().export().create("imgSlice", "Image");
        model.result().export("imgSlice").set("sourceobject", "pgImg");
        model.result().export("imgSlice").set("pngfilename",
                System.getProperty("user.dir") + java.io.File.separator + "figures"
                        + java.io.File.separator + "reverse_" + scheme + "_" + name + ".png");
        model.result().export("imgSlice").set("imagetype", "png");
        model.result().export("imgSlice").set("size", "manualweb");
        model.result().export("imgSlice").set("unit", "px");
        model.result().export("imgSlice").set("width", "1400");
        model.result().export("imgSlice").set("height", "700");
        model.result().export("imgSlice").set("resolution", "96");
        model.result().export("imgSlice").set("antialias", "on");
        model.result().export("imgSlice").set("zoomextents", "off");
        model.result().export("imgSlice").set("title", "on");
        model.result().export("imgSlice").set("legend", "on");
        model.result().export("imgSlice").set("logo", "off");
        model.result().export("imgSlice").set("options", "on");
        model.result().export("imgSlice").set("fontsize", "9");
        model.result().export("imgSlice").set("background", "color");
        model.result().export("imgSlice").set("customcolor", new double[]{1, 1, 1});
        model.result().export("imgSlice").run();
        model.result().export().remove("imgSlice");
        model.result().remove("pgImg");
        model.result().dataset().remove("cpImg");
    }

    private static void port(Model m, String tag, String sel, boolean excite, String name) {
        m.component("comp1").physics("ewfd").create(tag, "Port", 2);
        m.component("comp1").physics("ewfd").feature(tag).selection().named(sel);
        m.component("comp1").physics("ewfd").feature(tag).set("PortType", "Numeric");
        m.component("comp1").physics("ewfd").feature(tag).set("PortName", name);
        // explicit on BOTH ports: an excitation left at its default is not "off"
        m.component("comp1").physics("ewfd").feature(tag).set("PortExcitation", excite ? "on" : "off");
        if (excite) {
            m.component("comp1").physics("ewfd").feature(tag).set("Pin", 1);
        }
    }

    private static void blk(Model m, String tag, double x0, double y0, double z0,
                            double w, double h, double d) {
        m.component("comp1").geom("geom1").create(tag, "Block");
        m.component("comp1").geom("geom1").feature(tag).set("pos", new double[]{x0, y0, z0});
        m.component("comp1").geom("geom1").feature(tag).set("size", new double[]{w, h, d});
        m.component("comp1").geom("geom1").feature(tag).set("selresult", "on");
        m.component("comp1").geom("geom1").feature(tag).set("selresultshow", "all");
    }

    private static void matN(Model m, String tag, String sel, String n) {
        m.component("comp1").material().create(tag, "Common");
        m.component("comp1").material(tag).label("Mat " + tag);
        m.component("comp1").material(tag).selection().named(sel);
        m.component("comp1").material(tag).propertyGroup("def")
                .set("relpermittivity", new String[]{"(" + n + ")^2"});
        m.component("comp1").material(tag).propertyGroup("def")
                .set("relpermeability", new String[]{"1"});
        m.component("comp1").material(tag).propertyGroup("def")
                .set("electricconductivity", new String[]{"0"});
    }

    private static void boxSel(Model m, String tag, int dim, double x0, double x1,
                               double y0, double y1, double z0, double z1) {
        m.component("comp1").selection().create(tag, "Box");
        m.component("comp1").selection(tag).set("entitydim", dim);
        m.component("comp1").selection(tag).set("xmin", x0);
        m.component("comp1").selection(tag).set("xmax", x1);
        m.component("comp1").selection(tag).set("ymin", y0);
        m.component("comp1").selection(tag).set("ymax", y1);
        m.component("comp1").selection(tag).set("zmin", z0);
        m.component("comp1").selection(tag).set("zmax", z1);
        m.component("comp1").selection(tag).set("condition", "inside");
    }

    private static int[] union(Model m, String a, String b, String c) {
        int[] ea = m.component("comp1").selection(a).entities(2);
        int[] eb = m.component("comp1").selection(b).entities(2);
        int[] ec = m.component("comp1").selection(c).entities(2);
        int[] all = new int[ea.length + eb.length + ec.length];
        System.arraycopy(ea, 0, all, 0, ea.length);
        System.arraycopy(eb, 0, all, ea.length, eb.length);
        System.arraycopy(ec, 0, all, ea.length + eb.length, ec.length);
        return all;
    }

    private static void sizeDom(Model m, String tag, String sel, double h) {
        m.component("comp1").mesh("mesh1").create(tag, "Size");
        m.component("comp1").mesh("mesh1").feature(tag).selection().geom("geom1", 3);
        m.component("comp1").mesh("mesh1").feature(tag).selection().named(sel);
        m.component("comp1").mesh("mesh1").feature(tag).set("custom", "on");
        m.component("comp1").mesh("mesh1").feature(tag).set("hmaxactive", true);
        m.component("comp1").mesh("mesh1").feature(tag).set("hmax", String.valueOf(h));
    }

    private static void sizeBnd(Model m, String tag, String sel, double h) {
        m.component("comp1").mesh("mesh1").create(tag, "Size");
        m.component("comp1").mesh("mesh1").feature(tag).selection().geom("geom1", 2);
        m.component("comp1").mesh("mesh1").feature(tag).selection().named(sel);
        m.component("comp1").mesh("mesh1").feature(tag).set("custom", "on");
        m.component("comp1").mesh("mesh1").feature(tag).set("hmaxactive", true);
        m.component("comp1").mesh("mesh1").feature(tag).set("hmax", String.valueOf(h));
    }

    private static PrintWriter utf8(String name, boolean append) throws IOException {
        return new PrintWriter(new OutputStreamWriter(new FileOutputStream(name, append), "UTF-8"));
    }

    private static void write(String name, String line, boolean append) throws IOException {
        PrintWriter w = utf8(name, append);
        w.println(line);
        w.close();
    }
}
