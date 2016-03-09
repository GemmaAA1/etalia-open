define([
    'jquery',
    'app/ui/layout',
    'bootstrap'
], function($, layout) {

    function submitForm(event) {
        event.preventDefault();

        layout.setBusy();

        var $form = $(this);

        var submitXhr = $.ajax({
            url:  $form.attr('action'),
            type: $form.attr('method'),
            data: $form.serialize(),
            dataType: 'json'
        });

        submitXhr.done(function(data) {
            $('#id_errors').empty();

            // Redirect if key exists
            if (data.hasOwnProperty('redirect')) {
                window.location.href = data['redirect'];
            } else {
                checkUserInit();
            }
        });

        submitXhr.fail(function(data) {
            $('#id_errors').empty();
            $form.find('input').removeClass("alert alert-danger");

            var res = JSON.parse(data.responseText);
            $.each(res, function (key, value) {
                //console.log(key + ': ' + value);
                $form.find('input[name=' + key + ']')
                    .addClass("alert alert-danger");

                $form.find('#id_errors')
                    .prepend(
                        '<div class="error-message">' + '' +
                            '<span class="glyphicon glyphicon-exclamation-sign" aria-hidden="true"></span>' +
                            value +
                        '</div>'
                    );
            });

            layout.setAvailable();
        });

        return false;
    }

    function checkUserInit() {
        var statusInterval;

        layout.setBusy();

        statusInterval = setInterval(function() {
            $.getJSON('/user/user-update-step', function (data) {
                if (data.done) {
                    clearInterval(statusInterval);
                    if (data.hasOwnProperty('redirect')) {
                        window.location.href = data['redirect'];
                        return;
                    }
                    layout.setAvailable();
                } else {
                    // TODO improve response json format ...
                    layout.setBusy(
                        '<p><strong>' + data.messages[0] + '</strong></p>' +
                        '<p>' + data.messages[1] + '</p>'
                    );
                }
            });
        }, 1000);
    }

    $(function() {
        // Async form submission
        $('form[data-async]').on('submit', submitForm);
    });
});
