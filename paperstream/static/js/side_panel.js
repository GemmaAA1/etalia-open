/**
 * Created by nicolaspannetier on 10/19/15.
 */

$(document).ready(function() {

    $('.short').on('click', function () {
        if ($(this).closest('li').hasClass('active')) {  // unselect
            $(this).removeClass('active');
        } else { // select
            $(this).closest('li').addClass('active');
        }
        if ($(this).closest('li').hasClass('journal')){
            send_filter('partial', undefined);
        } else if ($(this).closest('li').hasClass('author')){
            send_filter(undefined, 'partial');
        }
    });
    $('.remove').on('click', function () {
        $(this).closest('li').removeClass('active');
        if ($(this).closest('li').hasClass('journal')){
            send_filter('partial', undefined);
        } else if ($(this).closest('li').hasClass('author')){
            send_filter(undefined, 'partial');
        }
    });
    $('.journal .only').on('click', function () {
        $('.journals').find('.active').each(function () {
            $(this).closest('li').removeClass('active')
        });
        $(this).closest('li').addClass('active');
        send_filter('partial', undefined);
    });
    $('.author .only').on('click', function () {
        $('.authors').find('.active').each(function () {
            $(this).closest('li').removeClass('active')
        });
        $(this).closest('li').addClass('active');
        send_filter(undefined, 'partial');
    });

    $('.journal, .author')
        .on('mouseenter', function () {
            $(this).find('.only').show();
            $(this).find('.short').show();
            $(this).find('.long').hide();
            $(this).find('.remove').show();
        })
        .on('mouseleave', function () {
            $(this).find('.only').hide();
            $(this).find('.short').hide();
            $(this).find('.long').show();
            $(this).find('.remove').hide();
        });

    $('.side-panel-title')
        .on('mouseenter', function (){
            $(this).find('.all').show();
        })
        .on('mouseleave', function () {
            $(this).find('.all').hide();
        });

    $('.journals .all').on('click', function () {
        console.log($(this).find('.journal'));
        $('.journals').find('.journal').each(function () {
            $(this).addClass('active')
        });
        send_filter('all', undefined);
    });

    $('.authors .all').on('click', function () {
        console.log($(this).find('.journal'));
        $('.authors').find('.author').each(function () {
            $(this).addClass('active')
        });
        send_filter(undefined, 'all');
    });

});

function send_filter(journal_flag, author_flag) {

    var json_data = {};
    json_data.action = 'filter';
    json_data.author_flag = author_flag;
    json_data.journal_flag = journal_flag;
    json_data.journals = [];
    json_data.authors = [];
    // build json object
    $('.journals').find('.active').each(function () {
        json_data['journals'].push($(this).data('journal-title'));
    });
    $('.authors').find('.active').each(function () {
        var author = {};
        author.first_name = $(this).data('first-name');
        author.last_name = $(this).data('last-name');
        json_data['authors'].push(author);
    });
    // ajax call and call back
    $.get(window.location.href, JSON.stringify(json_data), function (fragment) {
        $('.endless_data').html(fragment);
        $('.paper-list')
            .on('mouseenter', stamps_mouseenter)
            .on('mouseleave', stamps_mouseleave)
            .on('click', extendPaper).find('.no-toggling').click(function(event) {
                event.stopPropagation();
            });
        $('.tick').on('click', tick);
        $('.like').on('click', like);
        $('.add-to-library').on('click', add_to_lib);
        $('.trash').on('click', trash_paper);
    })
}