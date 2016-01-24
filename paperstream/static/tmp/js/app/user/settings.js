define(['jquery', 'app/ui/layout', 'bootstrap'], function($) {

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
                    console.log(json);
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
                    console.log(key + ': ' + value);
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

    $(function() {
        $('form[data-async]').on('submit', submitForm);
    });
});
