$(document).ready(function() {

    // bind envet / handlers for paper list
    bind_all();

    // send ajax with filter
    $('#most-recent').on('click', function (e) {
        e.preventDefault();
        send_filter(undefined, undefined, 'recent');
        $(this).closest('.dropdown-menu')
            .find('.sort-option')
            .removeClass('active');
        $(this).addClass('active');
    });

    $('#most-relevant').on('click', function (e) {
        e.preventDefault();
        send_filter(undefined, undefined, 'relevant');
        $(this).closest('.dropdown-menu')
            .find('.sort-option')
            .removeClass('active');
        $(this).addClass('active');
    });

    $('#most-trendy').on('click', function (e) {
        e.preventDefault();
        send_filter(undefined, undefined, 'trendy');
        $(this).closest('.dropdown-menu')
            .find('.sort-option')
            .removeClass('active');
        $(this).addClass('active');
    });

    $('#nothing-special').on('click', function (e) {
        e.preventDefault();
        send_filter(undefined, undefined, 'nothing');
        $(this).closest('.dropdown-menu')
            .find('.sort-option')
            .removeClass('active');
        $(this).addClass('active');
    });

});

function bind_all() {
    $('.paper-list')
        .on('mouseenter', stamps_mouseenter)
        .on('mouseleave', stamps_mouseleave)
        .on('click', extendPaper).find('.no-toggling').click(function(event) {
            event.stopPropagation();
        });

    // Short are liked paper that stay "active"
    $('.short').closest('.paper-list')
        .off('mouseleave')
        .off('mouseenter')
        .find('.stamps')
        .show();

    $('.feed .short, .trend .short')
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

    // Send trash paper from user library ajax call
    $('.restore').on('click', restore_paper);

    // tweet
    $('.tweet').on('click', share_tweet);

    // google plus
    $('.google-plus').on('click', share_google_plus);
}

// use for endless pagination callback
function unbind_all() {
    $('.paper-list')
        .off('mouseenter', stamps_mouseenter)
        .off('mouseleave', stamps_mouseleave)
        .off('click', extendPaper).find('.no-toggling').click(function(event) {
            event.stopPropagation();
        });

    // Send Like/tick ajax call
    $('.tick').off('click', tick);

    // Send Like/tick ajax call
    $('.like').off('click', like);

    // Send add paper to user library ajax call
    $('.add-to-library').off('click', add_to_lib);

    // Send trash paper from user library ajax call
    $('.trash').off('click', trash_paper);

    // tweet
    $('.tweet').off('click', share_tweet);

    // google plus
    $('.google-plus').off('click', share_google_plus);
}


function share_tweet () {
    var title = $(this).closest('.stamps-ext').data('title');
    var first_author = $(this).closest('.stamps-ext').data('first-compact');
    var url = $(this).closest('.stamps-ext').data('url');
    $.getJSON('http://api.bitly.com/v3/shorten?callback=?',
    {
        format: "json",
        apiKey: "R_8521db3604704f2ebaf044ec9be0ed6b",
        login: "nicolaspannetier",
        longUrl: url
    },
    function(response) {
        var text = title  + ' | ' + first_author + ' | ' + response.data.url + ' via @pubstreamio';
        //console.log(text);
        var geturl = 'https://twitter.com/intent/tweet?text=' + encodeURI(text);
        //console.log(geturl);
        //event.preventDefault();
        window.open(geturl, 'Tweet', 'width=600,height=400');
    });
}

function share_google_plus () {
    var title = $(this).closest('.stamps-ext').data('title');
    var first_author = $(this).closest('.stamps-ext').data('first-compact');
    var url = $(this).closest('.stamps-ext').data('url');
    var geturl = 'https://plus.google.com/share?url=' + encodeURI(url);
    event.preventDefault();
    window.open(geturl, 'Tweet', 'width=600,height=400');
}


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
        if (!$(this).children('.stamps').find('.like').hasClass('active')) {
            $(this).on('mouseenter', stamps_mouseenter)
                   .on('mouseleave', stamps_mouseleave);
        }
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
    var url = '/user/paper/ban';
    $.ajax({
        type: 'POST',
        url: url,
        data: {'pk': id,
               'source': window.location.pathname},
        success: function (json) {
            $.each(json, function (key, value) {
                if (key == 'is_banned') {
                    $tick.parents('.paper-list').slideUp(250);
                }
            });
        }
    });
}

function like (event) {
    var $like = $(this);
    $like.removeClass('like')
        .addClass('loading');
    var id = $(this).parents('.paper-list').attr('id');
    var url = '/user/paper/pin';
    $.ajax({
        type: 'POST',
        url: url,
        data: {'pk': id,
               'source': window.location.pathname},
        success: function (json) {
            $.each(json, function (key, value) {
                if (key == 'is_pinned') {
                    if (value) {
                        $like.removeClass('loading')
                            .addClass('like')
                            .addClass('active')
                            .parents('.paper-list')
                            .off('mouseleave', stamps_mouseleave)
                            .off('mouseenter', stamps_mouseenter)
                            .find('.stamps').show();
                    } else {
                        $like.removeClass('loading')
                            .addClass('like')
                            .removeClass('active')
                            .closest('.paper-list')
                            .on('mouseleave', stamps_mouseleave)
                            .on('mouseenter', stamps_mouseenter);
                    }
                } else {
                    $('#' + key).html(value);
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
                    }
                } else if (key == 'message') {
                    if (value) {
                        console.log(value);
                    }
                } else {
                    $('#' + key).html(value);
                }
            });
        }
    });
}

function restore_paper () {
    var $restore = $(this);
    $restore.removeClass('restore')
        .addClass('loading');
    var id = $(this).parents('.paper-list').attr('id');
    var url = '/user/paper/restore';
    $.ajax({
        type: 'POST',
        url: url,
        data: {pk: id},
        success: function (json) {
            $.each(json, function (key, value) {
                if (key == 'success') {
                    if (value) {
                        if(window.location.href.indexOf("library") > -1) {
                            $restore.parents('.paper-list').slideUp(250);
                        }
                    }
                } else if (key == 'message') {
                    if (value) {
                        console.log(value);
                    }
                } else {
                    $('#' + key).html(value);
                }
            });
        }
    });
}