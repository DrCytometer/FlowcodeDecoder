# FlowcodeDecoder

## A Shiny demultiplexer app for identifying FlowCode barcoded cells from flow cytometry data.

**Author: Orian Bricard**

FlowCodes are an adaptation of [ProCodes](https://www.cell.com/cell/fulltext/S0092-8674(18)31234-0?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS0092867418312340%3Fshowall%3Dtrue), a protein epitope barcoding system useful for CRISPR screening or other cell tracking experiments. 

The key difference between FlowCodes and ProCodes in that FlowCodes employ a histone backbone that allows for more stable and higher expression while preventing rejection in vivo. The app will work for both.

This Shiny app is intended to be an easy-to-use solution for debarcoding and identifying cells in a mixed population. In addition to identifying which FlowCode combination each cell has, the app permits the user to set gates and assess expression levels for any markers in the flow data set. Graphical outputs include summary plots of marker expression and pdfs of the thresholds for deciding barcode cutoffs. Data are exported in csv format both as summary data per sample and at the individual cell level for downstream processing. The thresholds (gates) can be re-used between analysis runs to facilitate processing in smaller batches or across experiments. You will need to export your cytometry data in CSV format with the biexponential (or arcsinh) scaling encoded (see exporting channel values in [FlowJo](https://docs.flowjo.com/flowjo/workspaces-and-samples/samples-and-file-types/ws-export/)). You will also need a csv file listing the barcodes and the meaningful label (e.g., CRISPR guide) that corresponds to each. See the example data for formatting.

For more information on FlowCodes, read the publication in [*Immunity*](https://www.cell.com/immunity/fulltext/S1074-7613(24)00277-2) or wait for the upcoming paper on this project. Researcher interested in using the FlowCode technology should contact Adrian Liston.
