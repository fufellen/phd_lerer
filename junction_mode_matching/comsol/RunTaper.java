import com.comsol.model.Model;
import com.comsol.model.util.ModelUtil;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.io.PrintWriter;
import java.util.Locale;

/**
 * Metal-width taper between the silicon feed (405 x 184 nm on the shared Sb2S3 buffer) and the
 * Sb2S3/Au phase shifter, with the TIP SHAPE as the parameter, driven from EITHER end.
 *
 * Geometry as in RunMatching.java (J_junction): half structure with a magnetic wall at x = 0,
 * the feed from z = -0.9 um to 0, the 450 nm silicon core from 0 to l_pcm, the gold as a
 * transition boundary condition on the footprint of a wedge body whose edge follows the profile
 * h(u), u = z/L, from a tip at z = 0 to the full half width at z = L.  Outside the footprint the
 * silicon sits straight on the buffer.
 *
 * Profiles (MT_SHAPE): linear, quad (parabola, sharp tip), sqrt (parabola, blunt tip), ellipse
 * (quarter ellipse, blunt tip), ellipsec (quarter ellipse, sharp tip), expo, gauss, rcos, klop,
 * fastwin, slowwin, step (control: no taper, metal starts abruptly at z = L), table (piecewise
 * linear profile from MT_TABLE = "u:h;u:h;...", h in units of the half width).
 * MT_SHAPE_L - taper length in um (default 0.5); MR_DIR = fwd | rev - which port is driven;
 * MR_MESH scales the mesh; MR_IMAGE exports slices; MR_LPCM overrides the phase-shifter length.
 * Every run writes the S-parameters of both ports and the port mode indices (taper_results.csv)
 * and the forward-power profile along z (taper_profiles.csv).
 */
public class RunTaper {

    private static final String STATUS = "taper_status.txt";
    private static PrintWriter out;
    private static PrintWriter prof;
    private static double meshScale = 1.0;
    private static double lFeed = 0.9;
    private static double lPcm = 3.0;
    private static String shape = "linear";
    private static double taperL = 0.5;
    private static double[][] table = null;
    private static int[] auFacesGlobal = new int[0];

    public static void main(String[] args) throws IOException {
        write(STATUS, "START taper", false);
        boolean fresh = !new java.io.File("taper_results.csv").exists();
        out = utf8("taper_results.csv", !fresh);
        if (fresh) {
            out.println("case,shape,taper_um,dir,n_pcm,mesh,l_feed_um,l_pcm_um,elements,"
                    + "neff1_re,neff1_im,neff2_re,neff2_im,"
                    + "S11_re,S11_im,S21_re,S21_im,S12_re,S12_im,S22_re,S22_im,p_feed_end,p_pcm_end");
        }
        out.flush();
        boolean freshP = !new java.io.File("taper_profiles.csv").exists();
        prof = utf8("taper_profiles.csv", !freshP);
        if (freshP) {
            prof.println("case,z_um,p_z");
        }
        prof.flush();

        meshScale = env("MR_MESH", 1.0);
        lFeed = env("MR_LFEED", 0.9);
        lPcm = env("MR_LPCM", 3.0);
        taperL = env("MT_SHAPE_L", 0.5);
        if (System.getenv("MT_SHAPE") != null) {
            shape = System.getenv("MT_SHAPE");
        }
        if (shape.equals("table")) {
            String t = System.getenv("MT_TABLE");
            if (t == null) {
                throw new IllegalArgumentException("MT_TABLE is required for shape=table");
            }
            String[] pts = t.split(";");
            table = new double[pts.length][2];
            for (int i = 0; i < pts.length; i++) {
                String[] uv = pts[i].split(":");
                table[i][0] = Double.parseDouble(uv[0].trim());
                table[i][1] = Double.parseDouble(uv[1].trim());
            }
        }
        String dir = System.getenv("MR_DIR");
        boolean forward = dir == null || dir.equals("fwd");
        double nPcm = env("MR_PCM", 2.712);
        boolean cryst = nPcm > 3.0;
        String shift2 = cryst ? "2.2772+1.50e-2*i" : "1.8877+2.46e-2*i";
        double shift1 = cryst ? env("MR_NFEED", 2.175) : 1.8877;
        run(forward, nPcm, shift1, shift2);
        out.close();
        prof.close();
        write(STATUS, "DONE", true);
    }

    private static double env(String name, double dflt) {
        String v = System.getenv(name);
        return v == null ? dflt : Double.parseDouble(v);
    }

    private static void run(boolean forward, double nPcm, double shift1, String shift2) throws IOException {
        String caseName = "sh_" + shape + "_L" + Math.round(taperL * 100) + "_" + (forward ? "fwd" : "rev")
                + (nPcm > 3.0 ? "_c" : "_a")
                + (meshScale != 1.0 ? "_m" + Math.round(meshScale * 100) : "")
                + (lPcm != 3.0 ? "_lp" + Math.round(lPcm * 100) : "");
        String tag = "MT_" + shape;
        write(STATUS, "CASE " + caseName, true);
        Model model = ModelUtil.create(tag);
        model.modelPath(System.getProperty("user.dir"));

        double wHalf = 1.0;
        double tPcm = 0.120, tLoad = 0.184356, wCore = 0.450, wFeed = 0.405;

        model.param().set("lambda0", "1.55[um]");
        model.param().set("f0", "c_const/lambda0");
        model.param().set("n_sio2", "1.444");
        model.param().set("n_si", "3.478");
        model.param().set("eps_au", "(0.6389-11.1748*i)^2");
        model.param().set("n_pcm", String.format(Locale.US, "%.4f", nPcm));
        model.param().set("t_au", "0.010[um]");

        model.component().create("comp1", true);
        model.component("comp1").geom().create("geom1", 3);
        model.component("comp1").geom("geom1").lengthUnit("um");
        blk(model, "bSub", 0, -tPcm - 0.70, -lFeed, wHalf, 0.70, lFeed + lPcm);
        blk(model, "bAir", 0, -tPcm, -lFeed, wHalf, 0.70 + tPcm, lFeed + lPcm);
        blk(model, "bPcm", 0, -tPcm, -lFeed, 0.400, tPcm, lPcm + lFeed);
        blk(model, "bWg", 0, 0.0, -lFeed, wFeed / 2.0, tLoad, lFeed);
        blk(model, "bCore0", 0, 0.0, 0.0, wCore / 2.0, tLoad, lPcm);
        curvedWedge(model, wCore / 2.0, tLoad, lPcm);
        model.component("comp1").geom("geom1").run();

        matN(model, "mAir", "geom1_bAir_dom", "1");
        matN(model, "mSub", "geom1_bSub_dom", "n_sio2");
        matN(model, "mPcm", "geom1_bPcm_dom", "n_pcm");
        matN(model, "mWg", "geom1_bWg_dom", "n_si");
        matN(model, "mCore0", "geom1_bCore0_dom", "n_si");
        matN(model, "mWedge", "geom1_hexAu_dom", "n_si");

        model.component("comp1").physics().create("ewfd", "ElectromagneticWavesFrequencyDomain", "geom1");
        boxSel(model, "selSym", 2, -1e-3, 1e-3, -1e4, 1e4, -1e4, 1e4);
        model.component("comp1").physics("ewfd").create("pmc1", "PerfectMagneticConductor", 2);
        model.component("comp1").physics("ewfd").feature("pmc1").selection().named("selSym");
        boxSel(model, "selOutX", 2, wHalf - 1e-3, wHalf + 1e-3, -1e4, 1e4, -1e4, 1e4);
        boxSel(model, "selOutYlo", 2, -1e4, 1e4, -tPcm - 0.70 - 1e-3, -tPcm - 0.70 + 1e-3, -1e4, 1e4);
        boxSel(model, "selOutYhi", 2, -1e4, 1e4, 0.70 - 1e-3, 0.70 + 1e-3, -1e4, 1e4);
        model.component("comp1").physics("ewfd").create("sbc1", "Scattering", 2);
        model.component("comp1").physics("ewfd").feature("sbc1").selection().set(
                union(model, "selOutX", "selOutYlo", "selOutYhi"));

        // the metal is the wedge body's own footprint on the plane y = 0
        boxSel(model, "selY0", 2, -1e-3, wCore / 2.0 + 1e-3, -1e-3, 1e-3, -1e-3, lPcm + 1e-3);
        model.component("comp1").selection().create("selWedge", "Intersection");
        model.component("comp1").selection("selWedge").set("entitydim", 2);
        model.component("comp1").selection("selWedge").set("input", new String[]{"geom1_hexAu_bnd", "selY0"});
        int[] auFaces = model.component("comp1").selection("selWedge").entities(2);
        model.component("comp1").selection().create("selAu", "Explicit");
        model.component("comp1").selection("selAu").geom(2);
        model.component("comp1").selection("selAu").set(auFaces);
        auFacesGlobal = auFaces;
        write(STATUS, "  wedge metal faces: " + auFaces.length, true);
        model.component("comp1").physics("ewfd").create("tbc1", "TransitionBoundaryCondition", 2);
        model.component("comp1").physics("ewfd").feature("tbc1").selection().named("selAu");
        model.component("comp1").physics("ewfd").feature("tbc1").set("d", "t_au");
        model.component("comp1").material().create("mAuSheet", "Common");
        model.component("comp1").material("mAuSheet").label("Gold sheet");
        model.component("comp1").material("mAuSheet").selection().geom("geom1", 2);
        model.component("comp1").material("mAuSheet").selection().named("selAu");
        model.component("comp1").material("mAuSheet").propertyGroup("def").set("relpermittivity", new String[]{"eps_au"});
        model.component("comp1").material("mAuSheet").propertyGroup("def").set("relpermeability", new String[]{"1"});
        model.component("comp1").material("mAuSheet").propertyGroup("def").set("electricconductivity", new String[]{"0"});

        boxSel(model, "selP1", 2, -1e4, 1e4, -1e4, 1e4, -lFeed - 1e-3, -lFeed + 1e-3);
        boxSel(model, "selP2", 2, -1e4, 1e4, -1e4, 1e4, lPcm - 1e-3, lPcm + 1e-3);
        port(model, "port1", "selP1", forward, "1");
        port(model, "port2", "selP2", !forward, "2");

        model.component("comp1").mesh().create("mesh1");
        model.component("comp1").mesh("mesh1").feature("size").set("custom", "on");
        model.component("comp1").mesh("mesh1").feature("size").set("hmax", String.valueOf(0.16 * meshScale));
        model.component("comp1").mesh("mesh1").feature("size").set("hmin", String.valueOf(0.006 * meshScale));
        model.component("comp1").mesh("mesh1").feature("size").set("hgrad", "1.5");
        sizeDom(model, "szWg", "geom1_bWg_dom", 0.070 * meshScale);
        sizeDom(model, "szPcm", "geom1_bPcm_dom", 0.060 * meshScale);
        sizeDom(model, "szCore0", "geom1_bCore0_dom", 0.060 * meshScale);
        sizeDom(model, "szWedge", "geom1_hexAu_dom", 0.055 * meshScale);
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

        double[] neff1 = global2(model, "ewfd.neff_1");
        double[] neff2 = global2(model, "ewfd.neff_2");
        double[] s11 = global2(model, "ewfd.S11");
        double[] s21 = global2(model, "ewfd.S21");
        double[] s12 = global2(model, "ewfd.S12");
        double[] s22 = global2(model, "ewfd.S22");
        write(STATUS, String.format(Locale.US,
                "  neff1=%.6f%+.3ei neff2=%.6f%+.3ei S11=%.4e%+.4ei S21=%.4e%+.4ei S12=%.4e%+.4ei S22=%.4e%+.4ei",
                neff1[0], neff1[1], neff2[0], neff2[1], s11[0], s11[1], s21[0], s21[1], s12[0], s12[1], s22[0], s22[1]), true);

        java.util.List<Double> zl = new java.util.ArrayList<Double>();
        zl.add(-lFeed + 0.1);
        zl.add(-0.30);
        zl.add(-0.10);
        for (double z = 0.25; z < lPcm - 0.15; z += 0.25) {
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
                "%s,%s,%.3f,%s,%.4f,%.2f,%.3f,%.3f,%d,%.7f,%.4e,%.7f,%.4e,"
                + "%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.9e,%.9e",
                caseName, shape, taperL, forward ? "fwd" : "rev", nPcm, meshScale, lFeed, lPcm, nel,
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
        write(STATUS, "CASE " + caseName + " OK", true);
        ModelUtil.remove(tag);
    }

    /** Half width of the metal along the taper as a fraction of its working value, u = z/L in [0, 1]. */
    private static double profile(double u) {
        if (u <= 0.0) {
            return 0.0;
        }
        if (u >= 1.0) {
            return 1.0;
        }
        if (shape.equals("linear")) {
            return u;
        }
        if (shape.equals("sqrt")) {
            return Math.sqrt(u);                              // parabola, blunt tip: width gained early
        }
        if (shape.equals("quad")) {
            return u * u;                                     // parabola, sharp tip: width gained late
        }
        if (shape.equals("quart")) {
            return Math.pow(u, 0.25);
        }
        if (shape.equals("cubic")) {
            return u * u * u;
        }
        if (shape.equals("ellipse")) {
            return Math.sqrt(1.0 - (1.0 - u) * (1.0 - u));    // quarter ellipse, blunt tip, flat at the end
        }
        if (shape.equals("ellipsec")) {
            return 1.0 - Math.sqrt(1.0 - u * u);              // quarter ellipse, sharp tip, steep at the end
        }
        if (shape.equals("rcos")) {
            return 0.5 * (1.0 - Math.cos(Math.PI * u));
        }
        if (shape.equals("expo")) {
            double a = 3.0;
            return (Math.exp(a * u) - 1.0) / (Math.exp(a) - 1.0);
        }
        if (shape.equals("gauss")) {
            double s = 0.35;
            double g0 = Math.exp(-1.0 / (2.0 * s * s));
            return (Math.exp(-(1.0 - u) * (1.0 - u) / (2.0 * s * s)) - g0) / (1.0 - g0);
        }
        if (shape.equals("klop")) {
            return klopfenstein(u);
        }
        if (shape.equals("fastwin")) {
            return piecewise(u, 0.45, 0.50);
        }
        if (shape.equals("slowwin")) {
            return piecewise(u, 0.25, 0.75);
        }
        if (shape.equals("table")) {
            for (int i = 0; i < table.length - 1; i++) {
                if (u >= table[i][0] && u <= table[i + 1][0]) {
                    double f = (u - table[i][0]) / (table[i + 1][0] - table[i][0]);
                    return table[i][1] + f * (table[i + 1][1] - table[i][1]);
                }
            }
            return u < table[0][0] ? table[0][1] : table[table.length - 1][1];
        }
        throw new IllegalArgumentException("unknown shape " + shape);
    }

    private static double piecewise(double u, double u0, double u1) {
        final double fLo = 210.0 / 450.0;
        final double fHi = 290.0 / 450.0;
        if (u <= u0) {
            return fLo * u / u0;
        }
        if (u <= u1) {
            return fLo + (fHi - fLo) * (u - u0) / (u1 - u0);
        }
        return fHi + (1.0 - fHi) * (u - u1) / (1.0 - u1);
    }

    private static double klopfenstein(double u) {
        final double a = 4.0;
        double raw = 0.5 + phi(2.0 * u - 1.0, a) / phi(1.0, a) * 0.5;
        double lo = 0.5 + phi(-1.0, a) / phi(1.0, a) * 0.5;
        double hi = 0.5 + phi(1.0, a) / phi(1.0, a) * 0.5;
        return (raw - lo) / (hi - lo);
    }

    private static double phi(double x, double a) {
        int n = 200;
        double s = 0.0;
        for (int i = 0; i < n; i++) {
            double y0 = x * i / n;
            double y1 = x * (i + 1) / n;
            s += 0.5 * (kern(y0, a) + kern(y1, a)) * (y1 - y0);
        }
        return s;
    }

    private static double kern(double y, double a) {
        double r = 1.0 - y * y;
        if (r <= 0.0) {
            return 0.5;
        }
        double t = a * Math.sqrt(r);
        return besselI1(t) / t;
    }

    private static double besselI1(double x) {
        double term = x / 2.0;
        double sum = term;
        for (int k = 1; k < 30; k++) {
            term *= (x * x / 4.0) / (k * (k + 1.0));
            sum += term;
        }
        return sum;
    }

    /** The wedge as an extruded polygon in the metal plane; "step" is the control with no taper. */
    private static void curvedWedge(Model model, double h, double tLoad, double lPcm) throws IOException {
        double[][] tab;
        if (shape.equals("step")) {
            tab = new double[][]{{taperL, 0.0}, {taperL, h}, {lPcm, h}, {lPcm, 0.0}};
        } else {
            // The polygon starts where the half width reaches H_MIN = 5 nm and is cut square there:
            // a sharp profile (quad, cubic, ellipsec) otherwise produces sub-nanometre slivers at the
            // tip, and the mesh of those slivers made one solve run for over an hour. A 10 nm wide
            // stub of metal is far below the 35 nm mesh on the sheet and changes nothing physical.
            final double hMin = 0.005;
            int n = Math.max(24, (int) Math.round(taperL / 0.020));
            // a TreeMap keyed by z sorts the points and merges duplicates; no anonymous comparator,
            // because comsolcompile emits only the top-level class (RunTaper$1 would be missing)
            java.util.TreeMap<Double, Double> curve = new java.util.TreeMap<Double, Double>();
            double uStart = -1.0;
            for (int i = 0; i <= 4000; i++) {
                double u = i / 4000.0;
                if (h * profile(u) >= hMin) {
                    uStart = u;
                    break;
                }
            }
            if (uStart < 0) {
                uStart = 1.0;
            }
            for (int i = 0; i <= n; i++) {
                double u = uStart + (1.0 - uStart) * i / n;
                curve.put(taperL * u, h * profile(u));
            }
            // blunt tips (sqrt, ellipse) rise steeply right after the cut: a few extra points there
            for (double du : new double[]{0.003, 0.008, 0.015}) {
                double u = uStart + du;
                if (u < 1.0) {
                    curve.put(taperL * u, h * profile(u));
                }
            }
            // square tip first, then the curve; drop points closer than 2 nm along z to their predecessor
            java.util.List<double[]> clean = new java.util.ArrayList<double[]>();
            clean.add(new double[]{taperL * uStart, 0.0});
            for (java.util.Map.Entry<Double, Double> e : curve.entrySet()) {
                double zc = e.getKey();
                double hc = e.getValue();
                if (clean.size() < 2 || zc - clean.get(clean.size() - 1)[0] > 0.002) {
                    clean.add(new double[]{zc, hc});
                }
            }
            tab = new double[clean.size() + 2][2];
            for (int i = 0; i < clean.size(); i++) {
                tab[i][0] = clean.get(i)[0];
                tab[i][1] = clean.get(i)[1];
            }
            tab[clean.size()][0] = lPcm;
            tab[clean.size()][1] = h;
            tab[clean.size() + 1][0] = lPcm;
            tab[clean.size() + 1][1] = 0.0;
        }
        model.component("comp1").geom("geom1").create("wpAu", "WorkPlane");
        model.component("comp1").geom("geom1").feature("wpAu").set("planetype", "quick");
        model.component("comp1").geom("geom1").feature("wpAu").set("quickplane", "zx");
        model.component("comp1").geom("geom1").feature("wpAu").set("quicky", 0.0);
        model.component("comp1").geom("geom1").feature("wpAu").geom().create("pol1", "Polygon");
        model.component("comp1").geom("geom1").feature("wpAu").geom().feature("pol1").set("source", "table");
        model.component("comp1").geom("geom1").feature("wpAu").geom().feature("pol1").set("type", "solid");
        model.component("comp1").geom("geom1").feature("wpAu").geom().feature("pol1").set("table", tab);
        model.component("comp1").geom("geom1").feature("wpAu").geom().run();
        model.component("comp1").geom("geom1").run("wpAu");
        model.component("comp1").geom("geom1").create("hexAu", "Extrude");
        model.component("comp1").geom("geom1").feature("hexAu").setIndex("distance", tLoad, 0);
        model.component("comp1").geom("geom1").feature("hexAu").selection("input").set("wpAu");
        model.component("comp1").geom("geom1").feature("hexAu").set("selresult", "on");
        model.component("comp1").geom("geom1").feature("hexAu").set("selresultshow", "all");
        write(STATUS, "  profile " + shape + ", points " + tab.length + ", taper " + taperL + " um", true);
    }

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

    private static void view(Model model, String scheme, String name, String plane, String offsetProp, double offset) {
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
        model.result().export("imgSlice").set("pngfilename", System.getProperty("user.dir") + java.io.File.separator
                + "figures" + java.io.File.separator + "taper_" + scheme + "_" + name + ".png");
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
        m.component("comp1").physics("ewfd").feature(tag).set("PortExcitation", excite ? "on" : "off");
        if (excite) {
            m.component("comp1").physics("ewfd").feature(tag).set("Pin", 1);
        }
    }

    private static void blk(Model m, String tag, double x0, double y0, double z0, double w, double h, double d) {
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
        m.component("comp1").material(tag).propertyGroup("def").set("relpermittivity", new String[]{"(" + n + ")^2"});
        m.component("comp1").material(tag).propertyGroup("def").set("relpermeability", new String[]{"1"});
        m.component("comp1").material(tag).propertyGroup("def").set("electricconductivity", new String[]{"0"});
    }

    private static void boxSel(Model m, String tag, int dim, double x0, double x1, double y0, double y1, double z0, double z1) {
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
