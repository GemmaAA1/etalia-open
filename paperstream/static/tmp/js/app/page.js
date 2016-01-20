define(['jquery', 'app/ui/layout', 'bootstrap'], function($, Layout) {

    $(function() {

        var layout = new Layout({debug: false});
        layout.init();

        $('#page-nav').affix();

        // FAQ
        var $faq = $('#faq');
        if ($faq.length) {
            $faq.on('click', 'a', function() {
                var target = $(this).attr('href');
                $faq.find('.collapse').not(target).collapse('hide');
                $(target).collapse('show');
            });
        }
    });
});
