$(document).ready(function() {

    $('.paper-list')
        .on('mouseenter', stamps_mouseenter)
        .on('mouseleave', stamps_mouseleave)
        .on('click', extendPaper).find('.no-toggling').click(function(event) {
            event.stopPropagation();
        });

    //
    $('.short').closest('.paper-list')
        .off('mouseleave')
        .off('mouseenter')
        .find('.stamps')
        .show();

    // Send Like/tick ajax call
    $('.tick').on('click', tick);

    // Send Like/tick ajax call
    $('.like').on('click', like);

    // Send add paper to user library ajax call
    $('.add-to-library').on('click', add_to_lib);

    // Send trash paper from user library ajax call
    $('.trash').on('click', trash_paper);

    // send ajax with filter
    $('#most-recent').on('click', function (e) {e.preventDefault(); send_filter(undefined, undefined, 'recent');});
    $('#most-relevant').on('click', function (e) {e.preventDefault(); send_filter(undefined, undefined, 'relevant');});
});

function stamps_mouseenter (){
    $(this).find('.stamps').show();
    $(this).find('.long').removeClass('long')
        .addClass('short');
}
function stamps_mouseleave () {
    $(this).find('.stamps').hide();
    $(this).find('.short').addClass('long')
        .removeClass('short');
}

function extendPaper(){
    if ($(this).children('.extended').is(':visible')) {
        $(this).children('.stamps').css('right', '25px');
        //$(this).children('.stamps').animate({"right":"22px"}, "fast");
        $(this).find('.more').removeClass('active');
        $(this).on('mouseenter', stamps_mouseenter)
               .on('mouseleave', stamps_mouseleave);

    } else {
        $(this).children('.stamps').css('right', '0px');
        //$(this).children('.stamps-div').animate({"right":"-10px"}, "fast");
        $(this).off('mouseenter mouseleave');
        $(this).find('.more').addClass('active');
    }
    $(this).children('.compact').toggle();
    $(this).children('.extended').toggle();
}

function tick () {
    var $tick = $(this);
    var id = $(this).parents('.paper-list').attr('id');
    var url = '/user/paper/tick';
    $.ajax({
        type: 'POST',
        url: url,
        data: {pk: id},
        success: function (json) {
            $.each(json, function (key, value) {
                if (key == 'is_ticked') {
                    $tick.parents('.paper-list').slideUp(250);
                }
            });
        }
    });
}

function like () {
    var $like = $(this);
    $like.removeClass('like')
        .addClass('loading');
    var id = $(this).parents('.paper-list').attr('id');
    var url = '/user/paper/like';
    $.ajax({
        type: 'POST',
        url: url,
        data: {pk: id},
        success: function (json) {
            $.each(json, function (key, value) {
                if (key == 'is_liked') {
                    if (value) {
                        $like.removeClass('loading')
                            .addClass('like')
                            .addClass('active')
                            .parents('.paper-list')
                            .off('mouseleave', stamps_mouseleave)
                            .off('mouseenter', stamps_mouseenter)
                            .find('.compact, .extended')
                            .addClass('bg-active');
                    } else {
                        $like.removeClass('loading')
                            .addClass('like')
                            .removeClass('active')
                            .closest('.paper-list')
                            .on('mouseleave', stamps_mouseleave)
                            .on('mouseenter', stamps_mouseenter)
                            .find('.compact, .extended')
                            .removeClass('bg-active');
                    }
                }
            });
        }
    });
}

function add_to_lib () {
    var $add = $(this);
    $add.removeClass('add-to-library')
        .addClass('loading');
    var id = $(this).parents('.paper-list').attr('id');
    var url = '/user/paper/add';
    $.ajax({
        type: 'POST',
        url: url,
        data: {pk: id},
        success: function (json) {
            $.each(json, function (key, value) {
                if (key == 'success') {
                    if (value) {  // success
                        $add.removeClass('loading')
                            .addClass('trash')
                            .attr('title', 'Move to trash')
                            .parents('.paper-list')
                            .addClass('bg-active')
                            .find('.like')
                            .addClass('active');
                        $add.on('click', trash_paper)
                            .off('click', add_to_lib);
                    } else {  // fail
                        $add.removeClass('loading')
                            .addClass('add-to-library');
                    }
                } else if (key == 'message') {
                    if (value) {
                        console.log(value);
                    }
                }
            });
        }
    });
}

function trash_paper () {
    var $trash = $(this);
    $trash.removeClass('trash')
        .addClass('loading');
    var id = $(this).parents('.paper-list').attr('id');
    var url = '/user/paper/trash';
    $.ajax({
        type: 'POST',
        url: url,
        data: {pk: id},
        success: function (json) {
            $.each(json, function (key, value) {
                if (key == 'success') {
                    if (value) {
                        if(window.location.href.indexOf("library") > -1) {
                            $trash.parents('.paper-list').slideUp(250);
                        } else {
                            $trash.removeClass('loading')
                                .addClass('add-to-library')
                                .attr('title', 'Add to your library');
                            $trash.on('click', add_to_lib)
                                .off('click', trash_paper);
                        }

                    } else {}
                } else if (key == 'message') {
                    if (value) {
                        console.log(value);
                    }
                }
            });
        }
    });
}