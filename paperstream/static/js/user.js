/**
 * Created by nicolaspannetier on 10/29/15.
 */

$(document).ready(function() {

    bind_profile_sections();
    bind_settings();

    // Update profile through modal ajax call
    $('.profile form[data-async]').on('submit', update_profile_forms);

    // Update settings through modal ajax call
    $('.settings form[data-async]').on('submit', update_settings_forms);

    // ajax call during signup
    $('.signup-forms form[data-async]').on('submit', update_signup_forms);


    // updates library
    $('#update-lib, #update-stream, #update-trend').on('click', function () {
        console.log($(this).attr('action'));
        $.ajax({
            url: $(this).attr('action'),
            success: function (resp) {
              $('#message').html(resp.message);
            },
            error: function (resp) {
                console.log('error');
                var res = JSON.parse(resp.responseText);
                $.each(res, function (key, value) {
                    $(this).siblings('.errors').html(value);
                });
            }
        });
    });

    $('.btn-file :file')
        .on('change', function() {
            var input = $(this),
                label = input.val().replace(/\\/g, '/').replace(/.*\//, '');
            input.parents('.avatar-upload-form').find('#file-name').html(label);
        });
});

function update_signup_forms(event) {
    var $form = $(this);
    $.ajax({
        type: $form.attr('method'),
        url: $form.attr('action'),
        data: $form.serialize(),

        success: function (json) {
            $('#id_errors').empty();
            //redirect if key exists
            if(json.hasOwnProperty('redirect')) {
                $(location).attr('href', json.redirect);
            }
        },

        error: function (resp) {
            console.log(resp.responseText);
            $('#id_errors').empty();
            var res = JSON.parse(resp.responseText);
            $('input').removeClass("alert alert-danger");
            $.each(res, function (key, value) {
                console.log(key + ': ' + value);
                $('input[name=' + key + ']').addClass("alert alert-danger");
                $('#id_errors').prepend('<div class="error-message"> <span class="glyphicon glyphicon-exclamation-sign" aria-hidden="true"></span>' + value + ' </div>');
            });
        }
    });
    event.preventDefault();
}

function update_profile_forms(event) {
    var $form = $(this);
    var $data = $($form.attr('data-target'));
    var $rootModal = $($form.attr('root-modal'));
    //console.log('form submitted!');
    $.ajax({
        type: $form.attr('method'),
        url: $form.attr('action'),
        data: $form.serialize(),

        success: function (json) {
            $('#id_errors').empty();
            // if all fields are empty switch to empty span else populate data
            // modal fields
            var all_fields_length = 0;
            $.each(json, function (key, value) {
                all_fields_length += value.length;
            });
            console.log(all_fields_length);
            if (all_fields_length == 0) {
                console.log($data.find('.has-value'));
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
            $rootModal.modal('hide');
            //redirect if key exists
            if(json.hasOwnProperty('redirect')) {
                $(location).attr('href', json.redirect);
            }
        },

        error: function (resp) {
            console.log(resp.responseText);
            $('#id_errors').empty();
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
}

function update_settings_forms(event) {
    var $form = $(this);
    var $data = $($form.attr('data-target'));
    var $rootModal = $($form.attr('root-modal'));
    //console.log('form submitted!');
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
                    $(location).attr('href', json.redirect);
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
                $('#id_errors').prepend('<div class="alert alert-danger" role="alert"> <span class="glyphicon glyphicon-exclamation-sign" aria-hidden="true"></span>' + value + ' </div>');
            });
        }
    });
    event.preventDefault();
}


function bind_profile_sections() {
    $('#affiliation-section, ' +
        '#names-section, ' +
        '#title-section, ' +
        '#position-section, ' +
        '#avatar-section')
        .on({
            'mouseenter': function (){
                $(this).find('a').css('color', 'rgba(127, 127, 127, 1)');
            },
            'mouseleave': function (){
                $(this).find('a').css('color', 'rgba(127, 127, 127, 0.4)');
            }
        });
}

function bind_settings() {
    $('#settings-section')
        .on({
            'mouseenter': function (){
                $(this).find('a').css('color', 'rgba(127, 127, 127, 1)');
            },
            'mouseleave': function (){
                $(this).find('a').css('color', 'rgba(127, 127, 127, 0.4)');
            }
        });
}

//function unbind_profile_sections() {
//    $('#affiliation-section, #names-section, #title-section, #position-section')
//        .off({
//            'mouseenter': function (){
//                $(this).find('a').css('color', 'rgba(127, 127, 127, 1)');
//            },
//            'mouseleave': function (){
//                $(this).find('a').css('color', 'rgba(127, 127, 127, 0.4)');
//            }
//        });
//}