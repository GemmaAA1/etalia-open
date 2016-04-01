define(['jquery', 'app/ui/layout', 'bootstrap', 'lib/d3'], function($) {
    // svg dimensions
    var w = 600;
    var h = 300;
    var padding = 30;

    var myData = 
        {articlesUploaded: 
            [{day: 10, articles: 20},
            {day: 20, articles: 14},
            {day: 30, articles: 20},
            {day: 40, articles: 21},
            {day: 50, articles: 15},
            {day: 60, articles: 22},
            {day: 70, articles: 9},
            {day: 80, articles: 6},
            {day: 90, articles: 23},
            {day: 100, articles: 7}]
        };

    function buildline(ds) {
        var xScale = d3.scale.linear() // or d3.time.scale()
            .domain([
                d3.min(ds.articlesUploaded, function(d){ return d.day; }),
                d3.max(ds.articlesUploaded, function(d){ return d.day; })
            ]) // input min/max
            .range([padding+5, w-padding]) // output px for svg
            .nice();
        var yScale = d3.scale.linear()
            .domain([0, d3.max(ds.articlesUploaded, function(d){ return d.articles; })])
            .range([h-padding, 10])
            .nice();

        var xAxisGen = d3.svg.axis().scale(xScale).orient("bottom");
        var yAxisGen = d3.svg.axis().scale(yScale).orient("left").ticks(4);

        var line = d3.svg.line()
            .x(function(d) { return xScale(d.day); })
            .y(function(d) { return yScale(d.articles); })
            .interpolate("basis"); //instead of linear
        
        var svg = d3.select(".data-viz").append("svg").attr({ width:w, height:h });
        
        var yAxis = svg.append("g").call(yAxisGen)
            .attr("class", "axis")
            .attr("transform", "translate(" + padding + ", 0)");

        var xAxis = svg.append("g").call(xAxisGen)
            .attr("class", "axis")
            .attr("transform", "translate(0," + (h-padding) + ")");

        var viz = svg.append("path")
            .attr({
                d: line(ds.articlesUploaded),
                "stroke": "grey",
                "stroke-width": 2,
                "fill": "none"
            });
    }

    // $.ajax({
        // method: 'GET',
        // url: '',
        // data: {
        //  format: 'json'
        // }})
        // .done(function(data) {
        //  console.log('data', data);
        // })
        // .fail(function(xrh, status, error) {
        //  console.log('failure', xrh, status, error);
        // });

    buildline(myData);

    
    function submitForm(e) {
        e.preventDefault();

        var $form = $(this),
            $data = $($form.data('target')),
            $modal = $($form.data('root'));

        //console.log('form submitted!');
        $.ajax({
            type: $form.attr('method'),
            url: $form.attr('action'),
            data: $form.serialize(),

            success: function (json) {
                $form.find('.form-errors').empty();

                // if all fields are empty switch to empty span else populate data
                // modal fields
                var all_fields_length = 0;
                $.each(json, function (key, value) {
                    all_fields_length += value.length;
                });
                //console.log(all_fields_length);

                if (all_fields_length == 0) {
                    //console.log($data.find('.has-value'));
                    $data.find('.has-value')
                        .addClass('hidden')
                        .siblings('.has-no-value')
                        .removeClass('hidden');
                } else {
                    $data.find('.has-value')
                        .removeClass('hidden')
                        .siblings('.has-no-value')
                        .addClass('hidden');

                    $.each(json, function (key, value) {
                        var $field = $('input[name=' + key + ']');
                        $field.val(value);
                        $field.removeClass("alert alert-danger");
                        $('#' + key).html(value);
                    });
                }
                $modal.modal('hide');

                //redirect if key exists
                if (json.hasOwnProperty('redirect')) {
                    window.location.href = json['redirect'];
                }
            },

            error: function (resp) {
                //console.log(resp.responseText);

                $form.find('.form-errors').empty();
                $('input').removeClass("alert alert-danger");

                var res = JSON.parse(resp.responseText);
                $.each(res, function (key, value) {
                    //console.log(key + ': ' + value);

                    $('input[name=' + key + ']').addClass("alert alert-danger");

                    $form.find('.form-errors').prepend(
                        '<div class="alert alert-danger" role="alert">' +
                            '<span class="glyphicon glyphicon-exclamation-sign" aria-hidden="true"></span>' +
                            value +
                        '</div>'
                    );
                });
            }
        });

        return false;
    }

    $(function() {
        $('form[data-async]').on('submit', submitForm);
    });
});
