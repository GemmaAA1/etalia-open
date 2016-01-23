$(document).ready(function() {

    $('.like-img').on('click', function () {
        $(this).toggleClass('active');
        send_filter(undefined, undefined);
    });

    $('.checkbox input').on('click', function () {
        send_filter();
    });

    $('#search-query').keypress(function (e) {
      if (e.which == 13) {
        send_filter();
        return false;
      }
    });

    //$('#search-query').on('input', function (e) {
    //    send_filter();
    //});


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

//function send_filter(journals_flag, authors_flag, sorting_flag) {
function send_filter() {

    var json_data = {};
    json_data.source = 'filter';
    json_data.like_flag = false;
    json_data.journals = [];
    json_data.authors = [];
    json_data.search_query = '';

    // build json object
    json_data.search_query = $('#search-query').val();

    $('.journals').find('input').each(function () {
        if ($(this).is(':checked')) {
            json_data['journals'].push($(this).data('journal-pk'));
            ga('send', 'event', 'Filter', 'select', $(this).data('journal-title'), parseInt($(this).data('journal-pk')));
        }
    });

    $('.authors').find('input').each(function () {
        if ($(this).is(':checked')) {
            json_data['authors'].push($(this).data('author-pk'));
            ga('send', 'event', 'Filter', 'select', $(this).data('author-fullname'), parseInt($(this).data('author-pk')));
        }
    });

    if ($('.like-img').hasClass('active')) {
        json_data.like_flag = true;
        ga('send', 'event', 'Pin', 'select');
    }

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