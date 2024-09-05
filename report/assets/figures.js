
var COLS = ["#1B98E0", "#df5353", "#ff9c33", "#f9d839", "#32a896"];
var COLS_AREA = ["#49afe9", "#e36969", "ffa84d"];
var COL_gray = "#666666";

// nreads vs chr length

var reads = alldata.dict_reads;

// var reads = {
//     'chrs': ['chr1',
//         'chr2',
//         'chr3',
//         'chr4',
//         'chr5',
//         'chr6',
//         'chr7',
//         'chr8',
//         'chr9',
//         'chr10',
//         'chr11',
//         'chr12',
//         'chr13',
//         'chr14',
//         'chr15',
//         'chr16',
//         'chr17',
//         'chr18',
//         'chr19',
//         'chr20',
//         'chr21',
//         'chr22',
//         'chrX',
//         'chrY',
//         'chrM',
//         'chrL'],
//     'lens': [248956422,
//         242193529,
//         198295559,
//         190214555,
//         181538259,
//         170805979,
//         159345973,
//         145138636,
//         138394717,
//         133797422,
//         135086622,
//         133275309,
//         114364328,
//         107043718,
//         101991189,
//         90338345,
//         83257441,
//         80373285,
//         58617616,
//         64444167,
//         46709983,
//         50818468,
//         156040895,
//         57227415,
//         16569,
//         48502],
//     'nreads': [66854410,
//         106747514,
//         61635942,
//         60943343,
//         52003515,
//         65196008,
//         43297414,
//         39560461,
//         32494057,
//         37695541,
//         34554069,
//         35915152,
//         30251220,
//         37098388,
//         20832687,
//         32931154,
//         23545670,
//         21700842,
//         10619470,
//         29545447,
//         12975906,
//         7782604,
//         21572724,
//         6667893,
//         5098632,
//         1937622]
// }

// sequence from 1
function seq_with(x){
    return Array.from(x, (val, i) => i + 1);
}

var updatemenus = [
    {
        buttons: [
            {
                args: [{'xaxis.type': 'linear', 'yaxis.type': 'linear'}],
                label: 'linear',
                method: 'relayout'
            },
            {
                args: [{'xaxis.type': 'log', 'yaxis.type': 'log'}],
                label: 'log',
                method: 'relayout'
            }
        ],
        font: {'size': 14},
        direction: 'left',
        pad: { 'r': 5, 't': 5 },
        showactive: true,
        type: 'buttons',
        x: 0.1,
        xanchor: 'left',
        y: 1.1,
        yanchor: 'top'
    }
]

var layout_reads_vs_lengths = {
    xaxis: {
        title: "chr/contig size (bp)",
        showgrid: false,
        showline: true,
        linewidth: 2,
        zeroline: false,
        ticks: "outside",
        tickwidth: 2,
        // type: 'log'
    },
    yaxis: {
        title: "reads mapped to each chr/contig",
        showgrid: false,
        showline: true,
        linewidth: 2,
        zeroline: false,
        ticks: "outside",
        tickwidth: 2,
        // type: 'log'
    },
    showlegend: false,
    margin: { b: 50, l: 90, t: 0, r: 0 },
    updatemenus: updatemenus
};

var config = { displayModeBar: false };

Plotly.newPlot(
    "plot-reads-vs-lengths",
    {
        data: [
            {
                x: reads.lens,
                y: reads.nreads,
                text: reads.chrs,
                name: "",
                mode: "markers",
                marker: { size: 5, color: COLS[0] },
            }
        ],
        layout: layout_reads_vs_lengths,
        config: config,
    }
);

