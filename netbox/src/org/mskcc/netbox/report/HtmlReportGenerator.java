package org.mskcc.netbox.report;

import edu.uci.ics.jung.algorithms.cluster.ClusterSet;
import edu.uci.ics.jung.algorithms.cluster.WeakComponentClusterer;
import edu.uci.ics.jung.graph.Graph;
import edu.uci.ics.jung.graph.Vertex;
import edu.uci.ics.jung.graph.decorators.StringLabeller;
import org.mskcc.netbox.algorithm.GeneConnector;
import org.mskcc.netbox.algorithm.LinkerGene;
import org.mskcc.netbox.algorithm.ModuleDetector;
import org.mskcc.netbox.genomic.GeneWithScore;
import org.mskcc.netbox.graph.GraphCreationException;
import org.mskcc.netbox.graph.GraphUtil;
import org.mskcc.netbox.graph.Module;
import org.mskcc.netbox.graph.NetworkStatsUtil;
import org.mskcc.netbox.graph.NewmanGirvanModuleDetector;
import org.mskcc.netbox.netcarto.NetCartoAnnealing;
import org.mskcc.netbox.util.Formatter;
import org.mskcc.netbox.util.GlobalConfig;
import org.mskcc.netbox.util.ProgressMonitor;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;

/**
 * HTML Report Generator.
 *
 * @author Ethan Cerami
 */
public final class HtmlReportGenerator {
    private static HtmlReportGenerator reportGenerator = new HtmlReportGenerator();
    private StringBuffer html = new StringBuffer();

    private HtmlReportGenerator() {
        a("<html>");
        a("<head>");
        a("<title>NetBox Report</title>");
        a("<style type=\"text/css\">");
        a("body    {\n"
                + "    margin: 20px;\n"
                + "    font-family: Verdana, Arial, sans-serif;\n"
                + "    color: black;\n"
                + "}\n");
        a("th {\n"
                + "    background: #AAAAAA;\n"
                + "    font-size:12px;\n"
                + "    text-align: left;\n"
                + "    color: white;\n"
                + "    padding:5px;\n"
                + "}\n"
                + "\n"
                + "td {\n"
                + "    font-family: Verdana, Arial, sans-serif;\n"
                + "    color: black;\n"
                + "    font-size:12px;\n"
                + "    margin:10px;\n"
                + "    margin-top:0px;\n"
                + "}");
        a(".green { background: #99FF99; padding:5px;}");
        a(".red { background: #FFCCCC; padding:5px;}");
        a("</style>");
        a("</head>");
        a("<body>");
        a("<h2>NetBox Report</h2>");
        GlobalConfig config = GlobalConfig.getInstance();

        String title = config.getTitle();
        if (title != null) {
            a("<h3>" + title + "</h3>");
        }

        a("<P>This report was auto-generated by NetBox on:  " + new Date());
        a("</p>");

        a("<h3>NetBox Parameters:</h3>");
        a("<table>");
        a("<tr><th>Parameter</th>");
        a("<th>Value</th>");
        a("</tr>");
        appendRow("Gene File", config.getGeneFile());
        appendRow("Mutation File", config.getMutationFile());
        appendRow("Copy Number Alteration File", config.getCnaFile());
        appendRow("Case File", config.getCaseFile());
        appendRow("Shortest Path Threshold", String.valueOf(config.getShortestPathThreshold()));
        if (config.getMutationFile() != null) {
            appendRow("Number of Random Trials", String.valueOf(config.getNumLocalTrials()));
            appendRow("Gene Frequency Threshold", String.valueOf(
                    config.getGeneFrequencyThreshold()));
            appendRow("Delta Threshold", String.valueOf(config.getDeltaThreshold()));
            appendRow("Module Frequency Threshold", String.valueOf(
                    config.getModuleFrequencyThreshold()));
            appendRow("Seed Radius Threshold", String.valueOf(config.getSeedRadiusThreshold()));
            if (config.getMRNAFile() != null) {
                appendRow("mRNA Z-Score Threshold", String.valueOf(
                        config.getMRNAZScoreThreshold()));
            }
            appendRow("Include Hemizygous Deletions and Single Copy Gains", String.valueOf(
                    config.includeLowLevelCnaChanges()));
        }
        int shortestPathThreshold = config.getShortestPathThreshold();
        if (shortestPathThreshold > 1) {
            appendRow("P-Value Linker Cut-Off", Formatter.getPValueFormat().format(
                    config.getPValueCutOff()));
        }
        String algo = config.getNetworkPartitionAlgorithm();
//        if (algo.equals(GlobalConfig.NG)) {
//            appendRow("Network Paritioning Algorithm", "Newman-Girvan");
//        } else {
//            appendRow("Network Paritioning Algorithm", "Simulated Annealing");
//        }
        a("</table>");
    }

    private void appendRow(String name, String value) {
        if (value != null) {
            a("<tr>");
            a("<td>" + name + "</td>");
            a("<td>" + value + "</td>");
            a("</tr>");
        }
    }

    private void appendRow(String name, File file) {
        if (file != null) {
            appendRow(name, file.getAbsolutePath());
        }
    }


    /**
     * Gets Singleton Instance.
     *
     * @return HTML Report Generator.
     */
    public static HtmlReportGenerator getInstance() {
        if (reportGenerator == null) {
            reportGenerator = new HtmlReportGenerator();
        }
        return reportGenerator;
    }

    /**
     * Appends Details about the Discovered Graph.
     *
     * @param geneList Gene List.
     * @param g        Graph Object.
     * @param geneConnector Gene Connector Object.
     */
    public void appendGraph(ArrayList<String> geneList, Graph g, GeneConnector geneConnector) {
        a("<h3>Network Discovered:  </h3>");
        a("<table>");
        a("<tr><th>Name</th>");
        a("<th>Value</th>");
        a("</tr>");
        appendRow("Number of input genes", String.valueOf(geneList.size()));
        appendRow("Number of vertices in graph", String.valueOf(g.numVertices()));
        appendRow("Number of edges in graph", String.valueOf(g.numEdges()));
        StringBuffer msg = new StringBuffer();
        if (GlobalConfig.getInstance().getShortestPathThreshold() == 1) {
            msg.append("At shortest path threshold of:  "
                    + GlobalConfig.getInstance().getShortestPathThreshold()
                    + ", I can "
                    + "connect " + geneConnector.getNumAlteredGenes() + " genes.");
        } else {
            msg.append("At shortest path threshold of:  "
                    + GlobalConfig.getInstance().getShortestPathThreshold()
                    + " and p-value cut-off of:  "
                    + Formatter.getPValueFormat().format(
                    GlobalConfig.getInstance().getPValueCutOff())
                    + ", I can "
                    + "connect " + geneConnector.getNumAlteredGenes() + " genes with "
                    + geneConnector.getNumLinkerGenes() + " linker genes.");
        }
        a("<tr><td colspan=2>" + msg.toString() + "</td></tr>");
        a("</table>");
    }

    /**
     * Appends Local Null Model Info.
     * @param s String.
     */
    public void appendLocalNullModel(String s) {
        a("<h3>Local Null Model:</h3>");
        a("<pre>");
        a(s);
        a("</pre>");
    }

    /**
     * Appends Global Null Model Info.
     * @param s String.
     */
    public void appendGlobalNullModel(String s) {
        a("<h3>Global Null Model:</h3>");
        a("<pre>");
        a(s);
        a("</pre>");
    }

    /**
     * Append Percent Module Detection Information.
     *
     * @param detector ModuleDetector detector.
     */
    public void appendPercentModuleDetection(ModuleDetector detector) {
        ArrayList<Module> moduleList = detector.getModules();

        if (moduleList.size() == 0) {
            a("<p>No Modules Detected.</p>");
        } else {
            a("<h3>Modules Detected:  </h3>");
            a("<table>");
            a("<tr>");
            a("<th>Module ID</th>");
            a("<th>Label</th>");
            a("<th>Frequency of Alteration</th>");
            a("<th>Genes</th>");
            a("<th># of Genes Explored</th>");
            a("<th>Unadjusted P-Value</th>");
            a("<th>FDR Adjusted P-Value</th>");
            a("</tr>");

            for (Module module : moduleList) {
                a("<tr>");
                a("<td>" + module.getModuleId() + "</td>");
                a("<td>" + module.getLabel() + "</td>");
                a("<td>" + Formatter.getShortDecimalFormat().format(module.getScore())
                        + "</td>");
                a("<td>");
                ArrayList<GeneWithScore> list = module.getGeneList();
                for (GeneWithScore gene : list) {
                    a(gene.getGene() + " ");
                }
                a("</td>");
                a("<td>" + (module.getGenesExplored().size() + 1) + "</td>");
                a("<td>" + outputPValue(module.getPValueUnAdjusted()) + "</td>");
                a("<td>" + outputPValue(module.getPValueFdrAdjusted()) + "</td>");
                a("</tr>");
            }
            a("</table>");

            a("<h3>Module Details:</h3>");
            for (Module module : moduleList) {
                ArrayList<GeneWithScore> list = module.getGeneList();
                if (list.size() > 1) {
                    a("<h4>Module ID:  " + module.getModuleId() + ", "
                            + module.getLabel() + "</h4>");
                    a("<table>");
                    a("<th>Gene Added in Order</th>");
                    a("<th>Gene Symbol</th>");
                    a("<th>Gene Set Frequency of Alteration</th>");
                    a("<th>Change</th>");
                    a("<tr>");

                    double currentScore = 0;
                    int order = 1;
                    for (GeneWithScore gene : list) {
                        a("<tr>");
                        a("<td>" + order + "</td>");
                        a("<td>" + gene.getGene() + "</td>");
                        a("<td>" + Formatter.getShortDecimalFormat().format(gene.getScore())
                                + "</td>");
                        double d = gene.getScore() - currentScore;
                        a("<td> +" + Formatter.getShortDecimalFormat().format(d) + "</td>");
                        a("</tr>");
                        currentScore = gene.getScore();
                        order++;
                    }
                    a("</table>");
                }
            }
        }
    }

    /**
     * Append Module Detection Information.
     *
     * @param detector   NewmanGirvanModuleDetector detector.
     * @param linkerList List of Linker Genes.
     */
    public void appendModuleDetection(NewmanGirvanModuleDetector detector,
                                      ArrayList<LinkerGene> linkerList) {
        Graph g = detector.getOptimalGraph();
        StringLabeller labeller = StringLabeller.getLabeller(g);
        WeakComponentClusterer wcSearch = new WeakComponentClusterer();
        ClusterSet moduleSet = wcSearch.extract(g);

        HashSet<String> linkerSet = new HashSet<String>();
        for (LinkerGene linker : linkerList) {
            linkerSet.add(linker.getGene());
        }

        if (moduleSet.size() > 0) {
            a("<h3>Modules Detected:  </h3>");
            a("<table>");
            a("<tr>");
            a("<th>Module ID</th>");
            a("<th>Number of Genes</th>");
            a("<th>Genes</th>");
            a("</tr>");

            for (int i = 0; i < moduleSet.size(); i++) {
                Graph gModule = moduleSet.getClusterAsNewSubGraph(i);
                a("<tr>");
                a("<td>" + i + "</td>");
                a("<td>" + gModule.getVertices().size() + "</td>");
                a("<td>");
                Iterator iterator = gModule.getVertices().iterator();
                while (iterator.hasNext()) {
                    Vertex vertex = (Vertex) iterator.next();
                    String label = GraphUtil.getVertexLabel(labeller, vertex);
                    if (linkerSet.contains(label)) {
                        label = label + "*";
                    }
                    a(label + " ");
                }
                a("</td>");
                a("</tr>");
            }
            if (linkerList.size() > 0) {
                a("<tr><td colspan=3>* Linker gene was not present in the original input list, "
                        + "but is "
                        + "significantly connected to members of the input list."
                        + "</td></tr>");
            }
            a("</table>");
        } else {
            a("<p>No modules detected.</p>");
        }
    }

    /**
     * Append Module Detection Information.
     *
     * @param detector   NewmanGirvanModuleDetector detector.
     * @param linkerList List of Linker Genes.
     */
    public void appendModuleDetection(NetCartoAnnealing detector,
                                      ArrayList<LinkerGene> linkerList) {
        HashSet<String> linkerSet = new HashSet<String>();
        for (LinkerGene linker : linkerList) {
            linkerSet.add(linker.getGene());
        }

        ArrayList<String> moduleList = detector.getFinalModuleList();
        HashMap<String, ArrayList<String>> moduleMap = detector.getGlobalModuleMap();
        if (moduleList.size() > 0) {
            a("<h3>Modules Detected:  </h3>");
            a("<table>");
            a("<tr>");
            a("<th>Module ID</th>");
            a("<th>Number of Genes</th>");
            a("<th>Genes</th>");
            a("</tr>");

            for (String moduleId : moduleList) {
                ArrayList<String> nodeList = moduleMap.get(moduleId);
                a("<tr>");
                a("<td>" + moduleId + "</td>");
                a("<td>" + nodeList.size() + "</td>");
                a("<td>");
                for (String node : nodeList) {
                    String label = node;
                    if (linkerSet.contains(node)) {
                        label = label + "*";
                    }
                    a(label + " ");
                }
                a("</td>");
                a("</tr>");
            }
            a("<tr><td colspan=3>* Linker gene was not present in the original input list, but is "
                    + "significantly connected to the network generated by the input list."
                    + "</td></tr>");
            a("</table>");
        } else {
            a("<p>No modules detected.</p>");
        }
    }


    /**
     * Append Linker Info.
     *
     * @param linkerList ArrayList of Linker Genes.
     * @throws GraphCreationException Graph Creation Error.
     */
    public void appendLinkers(ArrayList<LinkerGene> linkerList) throws GraphCreationException {
        if (linkerList.size() > 0) {
            a("<h3>Linker Gene Details:  Based on Global Network with:  "
                    + NetworkStatsUtil.getInstance().getNetworkStats().getNumGenes()
                    + " genes and "
                    + NetworkStatsUtil.getInstance().getNetworkStats().getNumEdges()
                    + " edges.</h3>");
            a("<table>");
            a("<tr>");
            a("<th>Gene Symbol</th>");
            a("<th>Local Degree</th>");
            a("<th>Global Degree</th>");
            a("<th>Unadjusted P-Value</th>");
            a("<th>FDR Adjusted P-Value</th>");
            a("<th>Status</th>");
            a("</tr>");

            for (LinkerGene linker : linkerList) {
                a("<tr>");
                a("<td>" + linker.getGene() + "</td>");
                a("<td>" + linker.getLocalDegree() + "</td>");
                a("<td>" + linker.getGlobalDegree() + "</td>");
                a("<td>" + outputPValue(linker.getUnadjustedPValue()) + "</td>");
                a("<td>" + outputPValue(linker.getFdrAdjustedPValue()) + "</td>");
                if (linker.getFdrAdjustedPValue()
                        > GlobalConfig.getInstance().getPValueCutOff()) {
                    a("<td class='red'>Pruned from Network</td>");
                } else {
                    a("<td class='green'>Included in Network</td>");
                }
                a("</tr>");
            }
        }
    }


    private String outputPValue(double p) {
        GlobalConfig globalConfig = GlobalConfig.getInstance();
        if (p > 0.0001) {
            NumberFormat formatter = new DecimalFormat("#0.0000");
            return formatter.format(p);
        } else if (p == 0) {
            double tempP = 1 / (double) globalConfig.getNumLocalTrials();
            NumberFormat formatter = new DecimalFormat("#0.0000");
            return "<" + formatter.format(tempP);
        } else {
            NumberFormat formatter = new DecimalFormat("0.#####E0");
            return formatter.format(p);
        }
    }

    /**
     * Finalize the HTML Report.
     *
     * @throws IOException IO Error.
     */
    public void finalizeReport() throws IOException {
        a("</body>");
        a("</html>");

        File file = new File("report.html");
        FileWriter writer = new FileWriter(file);

        ProgressMonitor.getInstance().setCurrentMessage("\nFinal HTML Report is available at:  "
                + file.getAbsolutePath());
        writer.write(html.toString());
        writer.close();
    }

    private void a(String m) {
        html.append(m + "\n");
    }
}