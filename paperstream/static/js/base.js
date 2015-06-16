jQuery(function ($) {

    $('form[data-async]').on('submit', function (event) {
        var $form = $(this);
        var $target = $($form.attr('data-target'));
        var $rootModal = $('div[root-modal]')
        console.log('form submitted!');
        $.ajax({
            type: $form.attr('method'),
            url: $form.attr('action'),
            data: $form.serialize(),

            success: function (json) {
                $('#id_errors').empty()
                $.each(json, function (key, value) {
                    var $field = $('input[name=' + key + ']');
                    $field.val(value);
                    $field.removeClass("alert alert-danger");
                });
                $rootModal.modal('hide');
                console.log('success');
            },

            error: function (resp) {
                console.log(resp.responseText);
                $('#id_errors').empty()
                var res = JSON.parse(resp.responseText);
                $('input').removeClass("alert alert-danger");
                $.each(res, function (key, value) {
                    console.log(key + ': ' + value);
                    $('input[name=' + key + ']').addClass("alert alert-danger");
                    $('#id_errors').prepend('<div class="alert alert-danger" role="alert"> <span class="glyphicon glyphicon-exclamation-sign" aria-hidden="true"></span>' + value + ' </div>');
                });
            }
        });

        event.preventDefault();
    });
});
