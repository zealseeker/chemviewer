function scatter_plot_run(data, cbs){
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
        .call(g => g.append("text")
            .attr("x", 4)
            .attr("y", margin.top)
            .attr("fill", "#000")
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
    let clip = svg.append("defs").append("SVG:clipPath")
        .attr("id", "clip")
        .append("SVG:rect")
        .attr("width", width )
        .attr("height", height )
        .attr("x", 0)
        .attr("y", 0);
    let scatter = svg.append('g')
        .attr('clip-path', 'url(#clip)'),
        xx = svg.append('g').call(xAxis),
        yy = svg.append('g').call(yAxis)
    scatter.selectAll(".dot")
        .data(data)
        .enter().append("circle")
        .attr("class", "dot")
        .attr("r", 3.5)
        .attr("transform", d => `translate(${x(d.PCA_0)},${y(d.PCA_1)})`)
        .style("fill", "#61a3a9")
        .style("opacity", 0.5)
        .on("click", clicked);
    let zoom = d3.zoom()
        .extent([[0, 0], [width, height]])
        .scaleExtent([0.5, 5])
        .on('zoom', zoomed)
        .on('end', zoomed)
    svg.append("rect")
        .attr("width", width)
        .attr("id", "selcmp")
        .attr("height", height)
        .style("fill", "none")
        .style("pointer-events", "all")
        .call(zoom);
    function clicked(event, d){
        d3.selectAll('.dot').attr('stroke', null)
        d3.select(this).attr("stroke", "red");
        if (cbs && cbs.clicked){
            cbs.clicked(d)
        }
    }
    let zoom_transform = {x: 0.0, y: 0.0};
    function zoomed({sourceEvent, transform}){
        var newX = transform.rescaleX(x);
        var newY = transform.rescaleY(y);
        // update axes with these new boundaries
        xx.call(d3.axisBottom(newX))
        yy.call(d3.axisLeft(newY))
        if (sourceEvent.type == "mousemove"){
            scatter.attr("transform", transform);
            let relative = {x: transform.x - zoom_transform.x,
                            y: transform.y - zoom_transform.y,
                            k: 1}
            scatter.attr('transform', `translate(${relative.x},${relative.y})`)
        }
        else{
            zoom_transform = transform
            scatter.attr('transform', null)
                .selectAll("circle")
                .attr('transform', d => `translate(${newX(d.PCA_0)},${newY(d.PCA_1)})`)
        }
     }
}
let scatter_plot = {
    run: scatter_plot_run,
    change_size(size){
        d3.select('div#scatter').selectAll('.dot').attr("r", size)
    },
    select_compounds(hide){
        let rect = d3.select('rect#selcmp')
        if (hide){
            rect.style('pointer-events', 'none')
        } else {
            rect.style('pointer-events', 'all')
        }
    }
}
export {scatter_plot}
