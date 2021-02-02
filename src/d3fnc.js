function scatter_plot_run(data){
    var margin = {top: 20, right: 20, bottom: 30, left: 40},
        width = 600 - margin.left - margin.right,
        height = 400 - margin.top - margin.bottom;
    var x = d3.scaleLinear()
        .domain(d3.extent(data, d=>d.PCA_0)).nice()
        .range([0, width]),
        y = d3.scaleLinear()
        .domain(d3.extent(data, d=>d.PCA_1)).nice()
        .range([height, 0])
    let yAxis = g => g
        .attr("transform", `translate(0,0)`)
        .call(d3.axisLeft(y))
        .call(g => g.select(".domain").remove())
        .call(g => g.select(".tick:last-of-type text").clone()
            .attr("x", 4)
            .attr("text-anchor", "start")
            .attr("font-weight", "bold")
            .text('PCA_1'))
    let xAxis = g => g
        .attr("transform", `translate(0,${height})`)
        .call(d3.axisBottom(x))
        .call(g => g.select(".domain").remove())
        .call(g => g.append("text")
            .attr("x", width - margin.right)
            .attr("y", -4)
            .attr("fill", "#000")
            .attr("font-weight", "bold")
            .attr("text-anchor", "end")
            .text('PCA_0'))
    let svg = d3.select('div#scatter').append('svg')
        .attr("width", width + margin.left + margin.right)
        .attr("height", height + margin.top + margin.bottom)
        .append("g")
        .attr("transform", "translate(" + margin.left + "," + margin.top + ")")
    
    svg.append("g")
        .call(xAxis);
    svg.append("g")
        .call(yAxis);

    svg.selectAll(".dot")
        .data(data)
        .enter().append("circle")
        .attr("class", "dot")
        .attr("r", 3.5)
        .attr("transform", d => `translate(${x(d.PCA_0)},${y(d.PCA_1)})`)
        .on("click", clicked);
    // svg.call(d3.zoom()
    //     .extent([[0, 0], [width, height]])
    //     .scaleExtent([0.5, 5])
    //     .on('zoom', zoomed)
    // )
    function clicked(event, d){
        d3.selectAll('.dot').attr('stroke', null)
        d3.select(this).attr("stroke", "red");
        console.log(d);
    }
    function zoomed({transform}){
        g.attr("transform", transform);
    }
}
let scatter_plot = {
    run: scatter_plot_run,
    change_size(size){
        d3.select('div#scatter').selectAll('.dot').attr("r", size)
    }
}
export {scatter_plot}
