/**
 * Created by nicolaspannetier on 10/19/15.
 */

$(document).ready(function() {

    $('.checkbox').on('click', function () {
        console.log('test');
        send_filter(undefined, undefined);
    });

    $('.journals .only').on('click', function (event) {
        $(this).closest('.side-panel-section').find('input').each(function () {
            $(this).prop('checked', false);
        });
        $(this).parents().siblings('input').prop('checked', true);
        send_filter('only', undefined);
        event.stopPropagation();
    });

    $('.authors .only').on('click', function (event) {
        $(this).closest('.side-panel-section').find('input').each(function () {
            $(this).prop('checked', false);
        });
        $(this).parents().siblings('input').prop('checked', true);
        send_filter(undefined, 'only');
        event.stopPropagation();
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

    //// show more block by block
    $('.show-more').on('click', function () {
         //get display block
        var last_vis_block = $(this).siblings('ul')
            .children('li:visible:last')
            .data('block');
        // hide show-more when reaching bottom
        var last_block = $(this).siblings('ul')
            .children('li:last')
            .data('block');
        var next = parseInt(last_vis_block) + 1;
        $(this).siblings('ul')
            .find("[data-block='" + next + "']")
            .show();

        if (last_block == next) {
            $(this).hide();
        }
         //show show-less
        $(this).siblings('.show-less').show();
    });

    // display show-less and bind onclick
    $('.show-less').on('click', function () {
        // get last vis block
        var last_vis_block = $(this).siblings('ul')
            .children('li:visible:last')
            .data('block');
        // hide block
        $(this).siblings('ul')
            .find("[data-block='" + last_vis_block + "']")
            .hide();
        // add show more if not visible
        $(this).siblings('.show-more').show();
        // remove show-less if reach top
        if (last_vis_block == "2") {
            $(this).hide();
        }
    });

});

function send_filter(journals_flag, authors_flag, sorting_flag) {

    var json_data = {};
    json_data.action = 'filter';
    json_data.journals_flag = journals_flag;
    json_data.authors_flag = authors_flag;
    json_data.sorting_flag = sorting_flag;
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
        bind_all();
        $.getScript( "https://d1bxh8uas1mnw7.cloudfront.net/assets/embed.js", function( data, textStatus, jqxhr ) {
                console.log('load performed');
        });
    })
}