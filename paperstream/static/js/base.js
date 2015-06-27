$(document).ready(function() {
    applyWhenElementExists('#user-lib-count-papers', '#loading-lib-block', update_lib_count, 1000);
});

function applyWhenElementExists(sel_to_up, sel_to_hide, myFunction, intervalTime) {

    var obj_up = jQuery(sel_to_up);
    var obj_hide = jQuery(sel_to_hide);
    var runningInterval = obj_up.data('interval');
    // if interval running, clear it
    if(runningInterval) {
        clearInterval(runningInterval);
        obj_up.removeData('interval');
    }
    // set interval
    var interval = setInterval(function() {
        if (obj_up.length > 0) {
            myFunction(obj_up, obj_hide);
        }
    }, intervalTime);
    // store interval in object for latter clearing in myFunction
    obj_up.data('interval', interval);
}

function update_lib_count(obj_up, obj_hide) {
    $.getJSON('/user/user-lib-count-papers/', function (json) {
        if (json.done) {
            var libInterval = obj_up.data('interval');
            clearInterval(libInterval);
            obj_up.removeData('interval');
            obj_hide.hide();
        }
        else {
            obj_up.html(json.count_papers);
        }
    });
}

jQuery(function ($) {
    // Profile update ajax call
    $('form[data-async]').on('submit', function (event) {
        var $form = $(this);
        var $target = $($form.attr('data-target'));
        var $rootModal = $($form.attr('root-modal'));
        //console.log('form submitted!');
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
                //redirect if key exists
                if(json.hasOwnProperty('redirect')) {
                    $(location).attr('href', json.redirect);
                }
                //console.log('success');
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


