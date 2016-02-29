define(['jquery', 'app/ui/layout', 'bootstrap'], function($) {

    $(function() {

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