/**
 * Created by nicolaspannetier on 11/4/15.
 */
$(document).ready(function() {
    // ajax call during signup
    $('.invite-form form[data-async]').on('submit', send_invite);
});

function send_invite(event) {
    var $form = $(this);
    var $data = $($form.attr('data-target'));
    var $rootModal = $($form.attr('root-modal'));
    $.ajax({
        type: $form.attr('method'),
        url: $form.attr('action'),
        data: $form.serialize(),
        success: function (json) {
            $.each(json, function (key, value) {
                // clean form
                $form.find("input[type=text], textarea").val("");
                // show success modal
                $rootModal.siblings('#invite-success')
                    .modal('show');
                // close modal
                $rootModal.modal('hide');
            })
        },
        error: function (resp) {
            console.log(resp.responseText);
            $('#id_errors').empty();
            var res = JSON.parse(resp.responseText);
            $.each(res, function (key, value) {
                $('#id_errors').prepend('<div class="alert alert-danger" role="alert"> <span class="glyphicon glyphicon-exclamation-sign" aria-hidden="true"></span>' + value + ' </div>');
            });
        }
    });
    event.preventDefault();
}
