
// flag to update locations of captions
var bin_img_is_updated = true;

function setWidth(id, dpiratio = 300 / 96) {
    let myImg = document.querySelector(id);
    let width = myImg.naturalWidth / dpiratio;
    widthnew = Math.min(930, Math.floor(width));
    $(id).width(widthnew + 'px');
}

function placeCaptions() {
    // $("#col-caption").empty();  // reset
    $(".caption").each(function (i, e) {
        $("#col-caption").append($(e));
        link = $(e).attr('link');
        pos = $(e).offset();
        pos.top = $(link).offset().top;
        $(e).offset(pos);
    });
}

function resetCaptions() {
    if (!bin_img_is_updated) return;
    $(".caption").each(function (i, e) {
        link = $(e).attr('link');
        pos = $(e).offset();
        pos.top = $(link).offset().top;
        $(e).offset(pos);
    });
    bin_img_is_updated = false;
}


// after page loading finished
$(function () {
    
    // add number marker to headers
    $("#main-content").children().filter('section').each(function (i, e) {
        h2 = $(e).find('h2').first();
        h2.text((i + 1) + '. ' + $(h2).text())
        $(e).find('h3').text(function (j, oritext) {
            return (i + 1) + '.' + (j + 1) + '. ' + oritext;
        });
    });

    $(".good").attr('title', 'primary');
    $(".neutral").attr('title', 'secondary');
    $(".bad").attr('title', 'last');

    // set initial size of figs
    setWidth("#img-depth-bin");
    setWidth("#img-meth-bin");
    if (alldata.mt_is_covered == 1) {
        setWidth("#cap-15", dpiratio = 300 / 110); // MT binning depth
    }
    if (alldata.plastid_is_covered == 1) {
        setWidth("#img-depth-bin-plastid", dpiratio = 300 / 110); // plastid binning depth
    }
    if (alldata.lambda_is_covered == 1) {
        setWidth("#cap-10", dpiratio=300 / 110); // lambda binning depth
    }

    // put the captions on the right panel
    placeCaptions();

    // methylation distribution
    function img_meth_dist() {
        dp = $("#range-dp").val();
        $("#label-1").text(dp);
        context = $("input[name='radio-1']:checked").val();
        file = "img/meth-dist-genome-" + context + "-dp" + dp + ".png";
        // console.log(file);
        $("#img-meth-dist").attr("src", file);
    }
    $("#range-dp").on("input propertychange", img_meth_dist);
    $("input[name='radio-1']").change(img_meth_dist);

    // chromosome wise meth
    function img_meth_chr() {
        dp = $("#range-2").val();
        $("#label-2").text(dp);
        context = $("input[name='radio-2']:checked").val();
        file = "img/meth-chr-" + context + "-double" + "-dp" + dp + ".png";
        // console.log(file);
        $("#img-meth-chr").attr("src", file);
    }
    $("#range-2").on("input propertychange", img_meth_chr);
    $("input[name='radio-2']").change(img_meth_chr);

    // binning  meth
    function img_meth_bin() {
        let dp = 1;
        let context = $("input[name='radio-3']:checked").val();
        let swicth = $("#switch-3").prop('checked');
        let strand = swicth ? 'single' : 'double';
        // file = "img/meth-bin-" + context + "-" + strand + "-dp" + dp + ".png";
        file = `img/meth-bin-${context}-${strand}-dp${dp}.png`;
        $("#img-meth-bin").attr("src", file);
        setWidth("#img-meth-bin");

        if (swicth) {
            $("#legend-meth-bin").show();
        } else {
            $("#legend-meth-bin").hide();
        }
        bin_img_is_updated = true;
    }
    $("#switch-3").on("change", img_meth_bin);
    $("input[name='radio-3']").on("change", img_meth_bin);

 
    // binning depth
    function img_depth_bin() {
        let swicth = $("#switch-4").prop('checked');
        let strand = swicth ? 'single' : 'double';
        file = `img/depth-bin-${strand}.png`;

        $("#img-depth-bin").attr("src", file);
        setWidth("#img-depth-bin");
        if (swicth) {
            $("#legend-depth-bin").show();
        } else {
            $("#legend-depth-bin").hide();
        }
        bin_img_is_updated = true;
    }    
    $("#switch-4").on("change", img_depth_bin);

    // meth of depth >= k
    function img_cul_meth_by_dp() {
        let context = $("input[name='radio-5']:checked").val();
        file = "img/meth-" + context + "-vs-dp-threshold.png";
        $("#img-cul-meth-by-dp").attr("src", file);
    }
    $("input[name='radio-5']").change(img_cul_meth_by_dp);

    // meth of depth = k
    function img_meth_at_dp() {
        let context = $("input[name='radio-meth-at-dp']:checked").val();
        file = `img/meth-${context}-at-dp-k.png`
        $("#img-meth-at-dp").attr("src", file);
    }
    $("input[name='radio-meth-at-dp']").change(img_meth_at_dp);

    // meth vs missing rate 
    function img_meth_vs_missing_rate() {
        let context = $("input[name='radio-meth-vs-missing-rate']:checked").val();
        file = `img/meth-${context}-vs-missing-rate.png`
        $("#img-meth-vs-missing-rate").attr("src", file);
    }
    $("input[name='radio-meth-vs-missing-rate']").change(img_meth_vs_missing_rate);

    // depth vs cytosine density
    function img_depth_vs_Cs_density() {
        let context = $("input[name='radio1-depth-vs-Cs-density']:checked").val();
        let strand = $("input[name='radio2-depth-vs-Cs-density']:checked").val();
        file = `img/depth-vs-${context}-density-of-${strand}-strand.png`
        $("#img-depth-vs-Cs-density").attr("src", file);
    }
    $("input[name='radio1-depth-vs-Cs-density']").change(img_depth_vs_Cs_density);
    $("input[name='radio2-depth-vs-Cs-density']").change(img_depth_vs_Cs_density);

    // meth vs missing rate 
    function img_meth_vs_missing_rate() {
        let context = $("input[name='radio-meth-vs-missing-rate']:checked").val();
        file = `img/meth-${context}-vs-missing-rate.png`
        $("#img-meth-vs-missing-rate").attr("src", file);
    }
    $("input[name='radio-meth-vs-missing-rate']").change(img_meth_vs_missing_rate);

    // // depth vs meth depth
    // function img_depth_vs_meth_depth() {
    //     let strand = $("input[name='radio-depth-vs-meth-depth']:checked").val();
    //     file = `img/depth-vs-meth-depth-of-${strand}-strand.png`
    //     $("#img-depth-vs-meth-depth").attr("src", file);
    // }
    // $("input[name='radio-depth-vs-meth-depth']").change(img_depth_vs_meth_depth);

    // depth AT vs CG
    function img_depth_AT_vs_CG() {
        let strand = $("input[name='radio-depth-AT-vs-CG']:checked").val();
        file = `img/depth-AT-vs-CG-of-${strand}-strand.png`
        $("#img-depth-AT-vs-CG").attr("src", file);
    }
    $("input[name='radio-depth-AT-vs-CG']").change(img_depth_AT_vs_CG);

    // meth vs cytosine density
    function img_meth_vs_Cs_density() {
        dp = $("#range-meth-vs-Cs-density").val();
        $("#label-meth-vs-Cs-density").text(dp);
        let context = $("input[name='radio1-meth-vs-Cs-density']:checked").val();
        let strand = $("input[name='radio2-meth-vs-Cs-density']:checked").val();
        file = `img/meth-vs-${context}-density-of-${strand}-strand-dp-ge${dp}.png`
        $("#img-meth-vs-Cs-density").attr("src", file);
    }
    $("input[name='radio1-meth-vs-Cs-density']").change(img_meth_vs_Cs_density);
    $("input[name='radio2-meth-vs-Cs-density']").change(img_meth_vs_Cs_density);
    $("#range-meth-vs-Cs-density").on("input propertychange", img_meth_vs_Cs_density);

    // coverage rate vs cytosine density
    function img_cov_rate_vs_Cs_density() {
        let dp = $("#range-covrate-vs-Cs-density").val();
        $("#label-covrate-vs-Cs-density").text(dp);
        let context = $("input[name='radio1-cov-rate-vs-Cs-density']:checked").val();
        let strand = $("input[name='radio2-cov-rate-vs-Cs-density']:checked").val();
        let file = `img/covrate-vs-${context}-density-of-${strand}-strand-dp-ge${dp}.png`;
        $("#img-cov-rate-vs-Cs-density").attr("src", file);
    }
    $("input[name='radio1-cov-rate-vs-Cs-density']").change(img_cov_rate_vs_Cs_density);
    $("input[name='radio2-cov-rate-vs-Cs-density']").change(img_cov_rate_vs_Cs_density);
    $("#range-covrate-vs-Cs-density").on("input propertychange", img_cov_rate_vs_Cs_density);

    // cytosine depth
    function img_cytosine_depth() {
        let context = $("input[name='radio-Cs-depth-dist']:checked").val();
        // file = `img/meth-${context}-at-dp-k.png`
        $("#img-Cs-depth-dist").attr("src", `img/genome-${context}-depth-distribution.png`);
        $("#img-Cs-depth-coverage").attr("src", `img/genome-${context}-coverage-vs-depth.png`);
    }
    $("input[name='radio-Cs-depth-dist']").change(img_cytosine_depth);

    // pangene coding
    function img_pangene_coding() {
        let context = $("input[name='radio-pangene-coding']:checked").val();
        $("#img-pan-gene-coding").attr("src", `img/pan-gene-meth-coding-${context}.png`);
    }
    $("input[name='radio-pangene-coding']").change(img_pangene_coding);

    // pangene lncRNA
    function img_pangene_lncRNA() {
        let context = $("input[name='radio-pangene-lncRNA']:checked").val();
        $("#img-pan-gene-lncRNA").attr("src", `img/pan-gene-meth-lncRNA-${context}.png`);
    }
    $("input[name='radio-pangene-lncRNA']").change(img_pangene_lncRNA);

    // pangene noncoding
    function img_pangene_noncoding() {
        let context = $("input[name='radio-pangene-noncoding']:checked").val();
        $("#img-pan-gene-noncoding").attr("src", `img/pan-gene-meth-other-noncoding-${context}.png`);
    }
    $("input[name='radio-pangene-noncoding']").change(img_pangene_noncoding);

    // update locations of captions, `img onload` fails
    document.querySelector('body').onwheel = resetCaptions;
    
});
