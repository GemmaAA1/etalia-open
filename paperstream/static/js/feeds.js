$(document).ready(function() {

    $('.paper-list')
        .on('mouseenter', stamps_mouseenter)
        .on('mouseleave', stamps_mouseleave)
        .on('click', extendPaper).find('.no-toggling').click(function(event) {
            event.stopPropagation();
        });

    // Send Like/tick ajax call
    $('.tick').on('click', tick);

    // Send Like/tick ajax call
    $('.like').on('click', like);

    // Send add paper to user library ajax call
    $('.add-to-library').on('click', add_to_lib);

    // Send trash paper from user library ajax call
    $('.trash').on('click', trash_paper);
});


function stamps_mouseenter (){
    $(this).find('.stamps').show();
    $(this).find('.short').show();
    $(this).find('.long').hide();
}
function stamps_mouseleave () {
    $(this).find('.stamps').hide();
    $(this).find('.short').hide();
    $(this).find('.long').show();
}

function extendPaper(){
    $(this).children('.compact').toggle();
    $(this).children('.extended').toggle();
    if ($(this).children('.extended').is(':visible')) {
        $(this).off('mouseenter mouseleave');
        $(this).find('.more').addClass('active');
    } else {
        $(this).find('.more').removeClass('active');
        $(this).on('mouseenter', stamps_mouseenter)
               .on('mouseleave', stamps_mouseleave);
        }
}

function tick (event) {
    var $tick = $(this);
    var id = $(this).parents('.paper-list').attr('id');
    var url = '/user/paper/tick';
    console.log(url);
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
    event.preventDefault();
}

function like (event) {
    var $like = $(this);
    var id = $(this).parents('.paper-list').attr('id');
    var url = '/user/paper/like';
    $.ajax({
        type: 'POST',
        url: url,
        data: {pk: id},
        success: function (json) {
            var tick = $(like).siblings('.tick');
            $.each(json, function (key, value) {
                if (key == 'is_liked') {
                    if (value) {
                        $like.addClass('active')
                            .parents('.paper-list')
                            .addClass('bg-active');
                    } else {
                        $like.removeClass('active')
                            .closest('.paper-list')
                            .removeClass('bg-active');
                    }
                }
            });
        }
    });
    event.preventDefault();
}

function add_to_lib (event) {
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
    //event.stopPropagation();
}

function trash_paper (event) {
    var $add = $(this);
    $add.removeClass('trash')
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
                    if (value) {
                        $add.removeClass('loading')
                            .addClass('add-to-library');
                        $add.on('click', add_to_lib)
                            .off('click', trash_paper);
                    } else {}
                } else if (key == 'message') {
                    if (value) {
                        console.log(value);
                    }
                }
            });
        }
    });
    event.stopPropagation();
}