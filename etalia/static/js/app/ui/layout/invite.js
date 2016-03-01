define(['jquery'], function ($) {

    var $inviteModal;

    function submitForm(event) {
        event.preventDefault();

        var $form = $inviteModal.find('.modal-form'),
            $result = $inviteModal.find('.modal-result');

        $.ajax({
            type: $form.attr('method'),
            url: $form.attr('action'),
            data: $form.serialize(),
            success: function (json) {
                // Clear errors
                $inviteModal.find('.form-errors').empty();

                // Clear controls
                $form.hide().find('input, textarea').val('');

                // Show result
                $result.show();
            },
            error: function (resp) {
                console.log(resp.responseText);

                var $errors = $inviteModal.find('.form-errors').empty();

                var res = JSON.parse(resp.responseText);
                $.each(res, function (key, value) {
                    $errors.prepend(
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

    $(function () {

        $inviteModal = $('#invite-modal');

        $('#open-invite').on('click', function() {
            $('#profile-dropdown').hide();
            $inviteModal.modal('show');
        });

        $inviteModal.find('form[data-async]').on('submit', submitForm);
        $inviteModal.find('.form-errors').empty();
        $inviteModal.on('hidden.bs.modal', function() {
            $inviteModal.find('.modal-form').show();
            $inviteModal.find('.modal-result').hide();
        });
    });
});
