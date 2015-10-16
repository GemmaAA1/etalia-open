$(document).ready(function() {

    toggleStamps();
    extendPaper();

    // Send Like/tick ajax call
    tick();

    // Send Like/tick ajax call
    like();

    //
    add_to_lib()
});


function toggleStamps(){
    $('.paper-list').on({
        mouseenter: function () {
            $(this).find('.stamps').show();
            $(this).find('.short').show();
            $(this).find('.long').hide();
        },
        mouseleave: function () {
            $(this).find('.stamps').hide();
            $(this).find('.short').hide();
            $(this).find('.long').show();
        }
    });
}

function extendPaper(){
    $('.paper-list').click( function () {
        console.log('toggle extend');
        $(this).children('.compact').toggle();
        $(this).children('.extended').toggle();
        if ($(this).children('.extended').is(':visible')) {
            console.log('toggle extend off');
            $(this).unbind('mouseenter mouseleave');
            $(this).find('.more').addClass('active');
        } else {
            console.log('toggle extend on');
            $(this).find('.more').removeClass('active');
            $(this).on({
                mouseenter: function () {
                    $(this).find('.stamps').show();
                    $(this).find('.stamps-lock').hide();
                },
                mouseleave: function () {
                    $(this).find('.stamps').hide();
                    $(this).find('.stamps-lock').show();
                }
            });
        }
    }).find('.open, .like, .tick, .go-to-pdf, .add-to-library' +
        '.comment, .email, .tweet, .google-plus').click(function(e) {
        e.stopPropagation();
    });
}

function tick () {
    $('.tick').on('click', function (event) {
        var tick = this;
        var id = $(this).parents('.paper-list').attr('id');
        var url = $(location).attr('protocol') + '//' +
            $(location).attr('host') + '/user/paper/tick';
        console.log(url);
        $.ajax({
            type: 'POST',
            url: url,
            data: {pk: id},
            success: function (json) {
                var like = $(tick).siblings('.like');
                $.each(json, function (key, value) {
                    if (key == 'is_liked') {
                        if (value) {
                            $(like).addClass('active');
                        } else {
                            $(like).removeClass('active');
                        }
                    } else if (key == 'is_ticked') {
                        $(tick).parents('.paper-list').slideUp(250);
                    }
                });
            }
        });
        event.preventDefault();
    });
}

function like () {
    $('.like').on('click', function (event) {
        var like = this;
        var id = $(this).parents('.paper-list').attr('id');
        var url = $(location).attr('protocol') + '//' +
            $(location).attr('host') + '/user/paper/like';
        $.ajax({
            type: 'POST',
            url: url,
            data: {pk: id},
            success: function (json) {
                var tick = $(like).siblings('.tick');
                $.each(json, function (key, value) {
                    if (key == 'is_liked') {
                        if (value) {
                            $(like).addClass('active');
                            $(like).parents('.paper-list').addClass('bg-active');
                        } else {
                            $(like).removeClass('active');
                            $(like).closest('.paper-list').removeClass('bg-active');
                        }
                    }
                });
            }
        });
        event.preventDefault();
    });
}

function add_to_lib () {
    $('.add-to-library').on('click', function (event) {
        console.log('add-to-lib');
        var add = this;
        var id = $(this).parents('.paper-list').attr('id');
        var url = $(location).attr('protocol') + '//' +
            $(location).attr('host') + '/user/paper/add';
        $.ajax({
            type: 'POST',
            url: url,
            data: {pk: id},
            success: function (json) {
                $.each(json, function (key, value) {
                    if (key == 'success') {
                        if (value) {
                            $(add).parents('.paper-list').addClass('bg-active');
                            $(add).parents('.paper-list').find('.like').addClass('active');
                            $(add).replaceWith('<span class="add-to-library active">Library</a>');
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
    });
}