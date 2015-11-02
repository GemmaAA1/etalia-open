$(document).ready(function() {

    applyWhenElementExists(
        $('#initializing-block'),
        $('#initializing-modal'),
        [$('#init-step'), $('#init-message')],
        '/user/user-init-step',
        update_init,
        1000);

    ////Modal message for lib syncing
    //applyWhenElementExists('#user-lib-count-matches', '#syncing-lib-block',
    //    '/user/user-lib-count-matches/', update_message, 1000);
    ////Modal message for feed updating
    //applyWhenElementExists('#user-feed-message', '#updating-feed-block',
    //    $(location).attr('href')+'user-feed-message', update_message, 2000);

    //// Profile update ajax call
    //$('form[data-async]').on('submit', function (event) {
    //    var $form = $(this);
    //    var $target = $($form.attr('data-target'));
    //    var $rootModal = $($form.attr('root-modal'));
    //    //console.log('form submitted!');
    //    $.ajax({
    //        type: $form.attr('method'),
    //        url: $form.attr('action'),
    //        data: $form.serialize(),
    //
    //        success: function (json) {
    //            $('#id_errors').empty();
    //            $.each(json, function (key, value) {
    //                var $field = $('input[name=' + key + ']');
    //                $field.val(value);
    //                $field.removeClass("alert alert-danger");
    //            });
    //            $rootModal.modal('hide');
    //            //redirect if key exists
    //            if(json.hasOwnProperty('redirect')) {
    //                $(location).attr('href', json.redirect);
    //            }
    //            //console.log('success');
    //        },
    //
    //        error: function (resp) {
    //            console.log(resp.responseText);
    //            $('#id_errors').empty();
    //            var res = JSON.parse(resp.responseText);
    //            $('input').removeClass("alert alert-danger");
    //            $.each(res, function (key, value) {
    //                console.log(key + ': ' + value);
    //                $('input[name=' + key + ']').addClass("alert alert-danger");
    //                $('#id_errors').prepend('<div class="alert alert-danger" role="alert"> <span class="glyphicon glyphicon-exclamation-sign" aria-hidden="true"></span>' + value + ' </div>');
    //            });
    //        }
    //    });
    //    event.preventDefault();
    //});


    $(".learn-more").click(function() {
       scrollToAnchor('learn-more');
    });
});


function applyWhenElementExists($block, $modal, $messages, url, myFunction, intervalTime) {
    // Control update of block / modal based on ajax calls

    var runningInterval = $block.data('interval');
    // if interval running, clear it
    if(runningInterval) {
        clearInterval(runningInterval);
        $block.removeData('interval');
    }
    // set interval
    var interval = setInterval(function() {
        if ($block.length > 0) {
            myFunction($block, $modal, $messages, url);
        }
    }, intervalTime);
    // store interval in object for latter clearing in myFunction
    $block.data('interval', interval);
}

function update_init($block,  $modal, $messages, url) {
    $.getJSON(url, function (json) {
        console.log(json);
        console.log($messages);
        if (json.done) {
            var libInterval = $block.data('interval');
            clearInterval(libInterval);
            $block.removeData('interval');
            $modal.hide()
                .siblings('#initializing-modal-backdrop')
                .hide();
        }
        else {
            $modal.show();
            $.each($messages, function (index, $sel) {
                console.log($sel);
                console.log(json.messages[index]);
                $sel.html(json.messages[index]);
            });
        }
    });
}


//function update_message(obj_up, obj_hide, url) {
//    $.getJSON(url, function (json) {
//        if (json.done) {
//            var libInterval = obj_up.data('interval');
//            clearInterval(libInterval);
//            obj_up.removeData('interval');
//            obj_hide.hide();
//            window.location.replace(json.url);
//        }
//        else {
//            obj_up.html(json.message);
//        }
//    });
//}

function scrollToAnchor(aid){
    var aTag = $("a[name='"+ aid +"']");
    $('html,body').animate({scrollTop: aTag.offset().top},'slow');
}