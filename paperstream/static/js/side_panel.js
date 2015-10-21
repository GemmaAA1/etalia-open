/**
 * Created by nicolaspannetier on 10/19/15.
 */

$(document).ready(function() {

    $('.checkbox').on('click', function () {
        send_filter(undefined, undefined);
    });

    $('.journals .only').on('click', function () {
        $(this).closest('.side-panel-section').find('input').each(function () {
            $(this).prop('checked', false);
        });
        $(this).parents().siblings('input').prop('checked', true);
        send_filter('only', undefined);
    });

    $('.authors .only').on('click', function () {
        $(this).closest('.side-panel-section').find('input').each(function () {
            $(this).prop('checked', false);
        });
        $(this).parents().siblings('input').prop('checked', true);
        send_filter(undefined, 'only');
    });

    $('.journals .all').on('click', function () {
        $(this).closest('.side-panel-section').find('input').each(function () {
            $(this).prop('checked', true);
        });
        send_filter('all', undefined);
    });

    $('.authors .all').on('click', function () {
        $(this).closest('.side-panel-section').find('input').each(function () {
            $(this).prop('checked', true);
        });
        send_filter(undefined, 'all');
    });

    $('.journal, .author')
        .on('mouseenter', function () {
            $(this).find('.short').show();
            $(this).find('.long').hide();
        })
        .on('mouseleave', function () {
            $(this).find('.short').hide();
            $(this).find('.long').show();
        });

    // display first block at load
    $('ul').find('[data-block="1"]').show();

    // show more block by block
    $('.show-more').on('click', function () {
        // get display block
        var last_vis_block = $(this).siblings('ul')
            .children('li:visible:last')
            .data('block');
        console.log(last_vis_block);
        // build next id to display
        var next = parseInt(last_vis_block) + 1;
        console.log(next);
        $(this).siblings('ul')
            .find("[data-block='" + next + "']")
            .show();

        // remove show-more if reach bottom
        var last_block = $(this).siblings('ul')
            .children('li:last')
            .data('block');
        if (last_block == next) {
            $(this).hide();
        }
    });


});

function send_filter(journals_flag, authors_flag) {

    var json_data = {};
    json_data.action = 'filter';
    json_data.journals_flag = journals_flag;
    json_data.authors_flag = authors_flag;
    json_data.journals = [];
    json_data.authors = [];

    // build json object
    $('.journals').find('input').each(function () {
        if ($(this).is(':checked')) {
            json_data['journals'].push([$(this).data('journal-title'), true]);
        }
        else {
            json_data['journals'].push([$(this).data('journal-title'), false]);
        }
    });

    $('.authors').find('input').each(function () {
        if ($(this).is(':checked')) {
            json_data['authors'].push([$(this).data('author-pk'), true]);
        }
        else {
            json_data['authors'].push([$(this).data('author-pk'), false]);
        }
    });

    console.log(JSON.stringify(json_data));

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