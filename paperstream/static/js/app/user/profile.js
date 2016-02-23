define(['jquery', 'app/ui/layout', 'bootstrap'], function($) {

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
