$(document).ready(function() {

    toggleStamps();
    extendPaper();

    // Send Like/tick ajax call
    tick();

    // Send Like/tick ajax call
    like();
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
        $(this).children('.compact').toggle();
        $(this).children('.extended').toggle();
        if ($(this).children('.extended').is(':visible')) {
            $(this).unbind('mouseenter mouseleave');
            $(this).find('.more').addClass('active');
        } else {
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
                var like2 = $(tick).parent().siblings('.compact').find('.like');
                var tick2 = $(tick).parent().siblings('.compact').find('.tick');
                console.log(like);
                console.log(like);
                $.each(json, function (key, value) {
                    if (key == 'is_liked') {
                        if (value) {
                            $(like).addClass('active');
                            $(like2).addClass('active');
                        } else {
                            $(like).removeClass('active');
                            $(like2).removeClass('active');
                        }
                    } else if (key == 'is_ticked') {
                        $(tick).parents('.paper-list').slideUp(250);
                        if (value) {
                            $(tick).addClass('active');
                            $(tick2).addClass('active');
                        } else {
                            $(tick).removeClass('active');
                            $(tick2).removeClass('active');
                        }
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
                var like2 = $(tick).parent().siblings('.compact').find('.like');
                var tick2 = $(tick).parent().siblings('.compact').find('.tick');
                $.each(json, function (key, value) {
                    if (key == 'is_liked') {
                        if (value) {
                            $(like).addClass('active');
                            $(like2).addClass('active');
                            console.log($(like).parents('.paper-list'));
                            $(like).parents('.paper-list').addClass('bg-active');
                        } else {
                            $(like).removeClass('active');
                            $(like2).removeClass('active');
                            $(like).closest('.paper-list').removeClass('bg-active');
                        }
                    } else if (key == 'is_ticked') {
                        if (value) {
                            $(tick).addClass('active');
                            $(tick2).addClass('active');
                        } else {
                            $(tick).removeClass('active');
                            $(tick2).removeClass('active');
                        }
                    }
                });
            }
        });
        event.preventDefault();
    });
}