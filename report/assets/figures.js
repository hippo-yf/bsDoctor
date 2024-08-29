
var COLS = ["#1B98E0", "#df5353", "#ff9c33", "#f9d839", "#32a896"];
var COLS_AREA = ["#49afe9", "#e36969", "ffa84d"];
var COL_gray = "#666666";

// freq on two strands

var data = alldata.cgkmer_2strands;

// sequence from 1
function seq_with(x){
    return Array.from(x, (val, i) => i + 1);
}

let minlim = Math.min(Math.min(...data.nW), Math.min(...data.nC));
let maxlim = Math.max(Math.max(...data.nW), Math.max(...data.nC));
// console.log(minlim, maxlim);

var layout_freq_on_two_strands = {
    // scattermode: "group",
    // title: "Grouped by Country",
    xaxis: {
        title: "relative freq on Watson strand",
        showgrid: false,
        showline: true,
        linewidth: 2,
        zeroline: false,
        ticks: "outside",
        tickwidth: 2,
    },
    yaxis: {
        title: "relative freq on Crick strand",
        showgrid: false,
        showline: true,
        linewidth: 2,
        zeroline: false,
        ticks: "outside",
        tickwidth: 2,
    },
    showlegend: false,
    margin: { b: 50, l: 50, t: 0, r: 0 },
};

var config = { displayModeBar: false };

Plotly.newPlot(
    "plot-freq-on-two-strands",
    {
        data: [
            {
                x: data.nW,
                y: data.nC,
                text: data.cgkmer,
                name: "",
                mode: "markers",
                marker: { size: 5, color: "#1B98E0" },
            },
            {
                x: [minlim, maxlim],
                y: [minlim, maxlim],
                mode: "lines",
                name: "",
                hoverinfo: "none",
                line: { color: COL_gray, lineWidth: 2, dash: "dash" },
            },
        ],
        layout: layout_freq_on_two_strands,
        config: config,
    }
);


// kmer meth on two strands

let minlim_medp1 = Math.min(Math.min(...data.meWdp1), Math.min(...data.meCdp1));
let maxlim_medp1 = Math.max(Math.max(...data.meWdp1), Math.max(...data.meCdp1));

var layout_me_on_two_strands = {
    xaxis: {
        title: "mean DNAme on Watson strand",
        showgrid: false,
        showline: true,
        linewidth: 2,
        zeroline: false,
        ticks: "outside",
        tickwidth: 2,
    },
    yaxis: {
        title: "mean DNAme on Crick strand",
        showgrid: false,
        showline: true,
        linewidth: 2,
        zeroline: false,
        ticks: "outside",
        tickwidth: 2,
    },
    showlegend: false,
    margin: { b: 50, l: 50, t: 0, r: 0 },
};

Plotly.newPlot(
    "plot-kmer-meth-on-two-strands",
    {
        data: [
            {
                x: data.meWdp1,
                y: data.meCdp1,
                text: data.cgkmer,
                name: "",
                mode: "markers",
                marker: { size: 5, color: "#1B98E0" },
            },
            {
                x: [minlim_medp1, maxlim_medp1],
                y: [minlim_medp1, maxlim_medp1],
                mode: "lines",
                name: "",
                hoverinfo: "none",
                line: { color: COL_gray, lineWidth: 2, dash: "dash" },
            },
        ],
        layout: layout_me_on_two_strands,
        config: config,
    }
);

// mean depth vs freq
let max_depth = Math.max(...data.mean_dp);
let min_depth = Math.min(...data.mean_dp);
var layout_depth_vs_freq = {
    // scattermode: "group",
    // title: "Grouped by Country",
    xaxis: {
        title: "relative freq on two strands",
        showgrid: false,
        showline: true,
        linewidth: 2,
        zeroline: false,
        ticks: "outside",
        tickwidth: 2,
    },
    yaxis: {
        title: "mean read depth",
        showgrid: false,
        showline: true,
        linewidth: 2,
        zeroline: false,
        ticks: "outside",
        tickwidth: 2,
    },
    showlegend: false,
    margin: { b: 50, l: 50, t: 0, r: 0 },
};

Plotly.newPlot(
    "plot-depth-vs-freq",
    {
        data: [
            {
                x: data.n,
                y: data.mean_dp,
                text: data.cgkmer,
                name: "",
                mode: "markers",
                marker: { size: 5, color: "#1B98E0" },
            },
            {
                x: [1, 1],
                y: [min_depth, max_depth],
                mode: "lines",
                name: "",
                hoverinfo: "none",
                line: { color: COL_gray, lineWidth: 2, dash: "dash" },
            },
        ],
        layout: layout_depth_vs_freq,
        config: config,
    }
);

// kmer meth vs freq
let max_me = Math.max(...data.medp1);
let min_me = Math.min(...data.medp1);
let layout_me_vs_freq = {
    xaxis: {
        title: "relative freq on two strands",
        showgrid: false,
        showline: true,
        linewidth: 2,
        zeroline: false,
        ticks: "outside",
        tickwidth: 2,
    },
    yaxis: {
        title: "mean DNAme level",
        showgrid: false,
        showline: true,
        linewidth: 2,
        zeroline: false,
        ticks: "outside",
        tickwidth: 2,
    },
    showlegend: false,
    margin: { b: 50, l: 50, t: 0, r: 0 },
};

Plotly.newPlot(
    "plot-kmer-meth-vs-freq",
    {
        data: [
            {
                x: data.n,
                y: data.medp1,
                text: data.cgkmer,
                name: "",
                mode: "markers",
                marker: { size: 5, color: "#1B98E0" },
            },
            {
                y: [min_me, max_me],
                x: [1, 1],
                mode: "lines",
                name: "",
                hoverinfo: "none",
                line: { color: COL_gray, lineWidth: 2, dash: "dash" },
            },
        ],
        layout: layout_me_vs_freq,
        config: config,
        // layout: { width: 600, height: 400 },
    }
);


function set_layout(xlabel, ylabel) {
    let layout = {
        xaxis: {
            title: xlabel,
            showgrid: false,
            showline: true,
            linewidth: 2,
            zeroline: false,
            ticks: "outside",
            tickwidth: 2,
        },
        yaxis: {
            title: ylabel,
            showgrid: false,
            showline: true,
            linewidth: 2,
            zeroline: false,
            ticks: "outside",
            tickwidth: 2,
        },
        showlegend: false,
        margin: { b: 50, l: 50, t: 0, r: 0 },
        hovermode: "closest",
    };
    return layout;
};

function hover_on(id) {
    return (data) => {
        let tn = data.points[0].curveNumber;
        let newcolor = "red";
        let update = {
            line: { color: newcolor, dash: "solid", width: 3 },
            opacity: 1,
        };
        Plotly.restyle(id, update, [tn]);
    };
};
function hover_off(id) {
    return (data) => {
        let tn = data.points[0].curveNumber;
        let oldcolor = "#1B98E0";
        let update = {
            line: { color: oldcolor, dash: "dot", width: 2 },
            opacity: 0.4,
        };
        Plotly.restyle(id, update, [tn]);
    };
};


// cgkmr cov vs depth
let data_plot_cov_vs_depth = [];
for (i = 0; i < data.cov.length; i++) {
    let trace = {
        // x: data.dp,
        x: seq_with(data.cov[i]),
        y: data.cov[i],
        name: "",
        text: data.cgkmer[i],
        mode: "lines",
        type: "scattergl",
        opacity: 0.4,
        line: { width: 2, color: "#1B98E0", dash: "dot" },
    };
    data_plot_cov_vs_depth.push(trace);
}

var myPlot = document.getElementById("plot-cov-vs-depth");

// var layout = {
//     // scattermode: "group",
//     // title: "Grouped by Country",
//     xaxis: {
//         title: "depth threshold",
//         showgrid: false,
//         showline: true,
//         linewidth: 2,
//         zeroline: false,
//         ticks: "outside",
//         tickwidth: 2,
//     },
//     yaxis: {
//         title: "coverage rate",
//         showgrid: false,
//         showline: true,
//         linewidth: 2,
//         zeroline: false,
//         ticks: "outside",
//         tickwidth: 2,
//     },
//     showlegend: false,
//     margin: { b: 50, l: 50, t: 0, r: 0 },
//     hovermode: "closest",
// };

var config = { displayModeBar: false };

Plotly.newPlot("plot-cov-vs-depth", data_plot_cov_vs_depth, set_layout(xlabel = "depth threshold", ylabel = "coverage rate"), config);
myPlot.on("plotly_hover", hover_on("plot-cov-vs-depth"));
myPlot.on("plotly_unhover", hover_off("plot-cov-vs-depth"));


// cgkmr me vs depth
let data_plot_me_vs_depth = [];
console.log(data.me.length);

for (i = 0; i < data.me.length; i++) {
    let trace = {
        // x: data.dp,
        x: seq_with(data.me[i]),
        y: data.me[i],
        name: "",
        text: data.cgkmer[i],
        mode: "lines",
        type: "scattergl",
        opacity: 0.4,
        line: { width: 2, color: "#1B98E0", dash: "dot" },
    };
    data_plot_me_vs_depth.push(trace);
}

var myPlot = document.getElementById("plot-me-vs-depth");

Plotly.newPlot("plot-me-vs-depth", data_plot_me_vs_depth, set_layout(xlabel = "depth threshold", ylabel = "mean DNAme level"), config);
myPlot.on("plotly_hover", hover_on("plot-me-vs-depth"));
myPlot.on("plotly_unhover", hover_off("plot-me-vs-depth"));
