define(['jquery', 'app/ui/layout', 'jquery-ui', 'bootstrap'], function($, layout) {

    function submitForm(e) {
        var $form = $(this);
        var $data = $($form.attr('data-target'));
        var $rootModal = $($form.attr('root-modal'));

        console.log('form submitted!');
        console.log($data);

        $.ajax({
            type: $form.attr('method'),
            url: $form.attr('action'),
            data: $form.serialize(),
            success: function (json) {
                $.each(json, function (key, value) {
                    //console.log(json);
                    var $field = $('input[name=' + key + ']');
                    $field.val(value);
                    $field.removeClass("alert alert-danger");
                    $('#' + key).html(value);
                    $rootModal.modal('hide');
                    //redirect if key exists
                    if (json.hasOwnProperty('redirect')) {
                        $(location).attr('href', json['redirect']);
                    }
                })
            },
            error: function (resp) {
                console.log(resp.responseText);
                $('#id_errors').empty();
                var res = JSON.parse(resp.responseText);
                $('input').removeClass("alert alert-danger");
                $.each(res, function (key, value) {
                    //console.log(key + ': ' + value);
                    $('input[name=' + key + ']').addClass("alert alert-danger");
                    $('#id_errors').prepend(
                        '<div class="alert alert-danger" role="alert">' +
                            '<span class="glyphicon glyphicon-exclamation-sign" aria-hidden="true"></span>' + value +
                        '</div>'
                    );
                });
            }
        });

        e.preventDefault();
        return false;
    }

    function initSliders() {
        $('#id_stream_vector_weight, ' +
            '#id_stream_author_weight, ' +
            '#id_stream_journal_weight, ' +
            '#id_trend_doc_weight, ' +
            '#id_trend_altmetric_weight')
            .hide()
            .wrap('<div class="slider"></div>');

        $('.slider').each(function(i, slider) {
            var $slider = $(slider),
                $input = $slider.find('input');

            $slider.slider({
                range: false,
                animate: false,
                min: parseFloat($input.data('slider-min'))*100,
                max: parseFloat($input.data('slider-max'))*100,
                step: parseFloat($input.data('slider-step'))*100,
                value: parseFloat($input.val())*100,
                slide: function( event, ui ) {
                    var value = ui.value/100;

                    $slider.find('span')
                        .tooltip('destroy')
                        .tooltip({
                            title: value,
                            animation: false,
                            container: $slider.find('span')
                        })
                        .tooltip('show');

                    $input.val(value);
                },
                create: function() {
                    $slider.find('span')
                        .tooltip({
                            title: $input.val(),
                            animation: false,
                            container: $slider.find('span')
                        });
                }
            });
        });
    }

    $(function() {
        // ASync form submission
        $('form[data-async]').on('submit', submitForm);

        // Jquery UI Sliders
        initSliders();


        $('#update-stream').on('click', function() {
            layout.setBusy('<p>Test message</p>');
        });
    });
});
