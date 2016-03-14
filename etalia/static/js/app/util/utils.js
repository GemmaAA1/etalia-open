define(['jquery'], function($) {

    // CSRF helper functions taken directly from Django docs
    function getCookie(name) {
        var cookieValue = null;
        if (document.cookie && document.cookie != '') {
            var cookies = document.cookie.split(';');
            for (var i = 0; i < cookies.length; i++) {
                var cookie = jQuery.trim(cookies[i]);
                // Does this cookie string begin with the name we want?
                if (cookie.substring(0, name.length + 1) == (name + '=')) {
                    cookieValue = decodeURIComponent(cookie.substring(name.length + 1));
                    break;
                }
            }
        }
        return cookieValue;
    }
    var csrftoken = getCookie('csrftoken');
    function csrfSafeMethod(method) {
        // these HTTP methods do not require CSRF protection
        return (/^(GET|HEAD|OPTIONS|TRACE)$/.test(method));
    }
    $.ajaxSetup({
        beforeSend: function (xhr, settings) {
            if (!csrfSafeMethod(settings.type) && !this.crossDomain) {
                xhr.setRequestHeader("X-CSRFToken", csrftoken);
            }
        }
    });



    return {
        toggleClass: function ($element, cssClass, state) {
            if (typeof state != 'undefined') {
                if (state) {
                    $element.addClass(cssClass);
                } else {
                    $element.removeClass(cssClass);
                }
                return state;
            }
            if ($element.hasClass(cssClass)) {
                $element.removeClass(cssClass);
                return false;
            }
            $element.addClass(cssClass);
            return true;
        },
        getParameterByName: function (name, url) {
            if (!url) url = window.location.href;
            name = name.replace(/[\[\]]/g, "\\$&");
            var regex = new RegExp("[?&]" + name + "(=([^&#]*)|&|#|$)"),
                results = regex.exec(url);
            if (!results) return null;
            if (!results[2]) return '';
            return decodeURIComponent(results[2].replace(/\+/g, " "));
        },
        bindLoadingButtons: function($element) {
            $element.on('click', '.loading-button', function(e) {
                var $button = $(e.target).closest('.loading-button');

                if ($button.data('busy')) {
                    e.stopPropagation();
                    e.stopImmediatePropagation();
                    return;
                }

                $button.data('busy', true);

                $button.find('.eai')
                    .removeAttr('class')
                    .addClass('eai eai-loading');
            });
        },
        restoreLoadingButton: function($button, iconClass) {
            $button.data('busy', false);

            $button.find('.eai')
                .removeAttr('class')
                .addClass('eai ' + iconClass);
        },
        popup: function(url, name, width, height) {
            name = name || 'etalia-popup';
            width = width || 520;
            height = height || 460;

            var left = Math.round((window.innerWidth/ 2) - (width / 2)),
                top = Math.round((window.innerHeight / 2) - (height / 2)),
                params = "menubar=no,toolbar=no,resizable=yes,scrollbars=yes," +
                         "width=" + width + ",height=" + height + "," +
                         "top=" + top + ",left=" + left;

            return window.open(url, name, params);
        }
    }
});